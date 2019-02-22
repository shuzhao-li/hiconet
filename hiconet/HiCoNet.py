# Copyright 2018 Shuzhao Li. All Rights Reserved.
#
# Licensed under the BSD 3-Clause License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Hieracrchical Community Network
for integration of multiple data types collected from a common group of subjects.

Terminology:
    Study: Same as ImmPort "Study"

    Project: a collection of data of one or more types. For multiple data types, common samples/subjects are expected.
        This is the unit HiCoNet works on - HiCoNet integrates DataMatrices within a DataSet
        A DataSet should have at least one Society of data


    Society: one data type
    at least one of each DataMatrix, FeatureAnnotation, ObservationAnnotation

    DataMatrix: a data matrix of continuous values that represent a biological state or concentration, of the same data type.
        This can include different time points or treatments.
        This is the unit community detection is based on.

    ObservationAnnotation: meta data on samples. This may include TimePoints and Treatments, often in biosample table from ImmPort DB
    FeatureAnnotation: meta data on features

    Key annotation variables: time point and treatment.
    TimePoint:
    Treatment:


    Graph: as in loom (loompy.org).
    Community: a group of features within a society that share a similar pattern.
    e.g. one way to identify communities is to do Louvain graph-clustering method on the graph.


    Reusing
    -------
    The data structure is aligned with anndata (https://github.com/theislab/anndata)
    For loom format, one loom file = a Society; as meta data can differ for different data types.
    Same goes for anndata


Matching samples/subjects:
    DataSets are rarely complete. Sample matching is done on the fly when needed.
The biosample table typically found in ImmPort:
    BIOSAMPLE_ACCESSION    DESCRIPTION    NAME    PLANNED_VISIT_ACCESSION    STUDY_ACCESSION    STUDY_TIME_COLLECTED    STUDY_TIME_COLLECTED_UNIT    STUDY_TIME_T0_EVENT    STUDY_TIME_T0_EVENT_SPECIFY    SUBJECT_ACCESSION    SUBTYPE    TYPE    WORKSPACE_ID

default integration schema:
    btw all data types and all time points
    users can specifiy the schema too.
    
    If data are normalized to baseline values, a new society is derived from old society.

        ? Requiring _Minimal_Sample_Number,
        ? Check of data quality and imputation will be added later.

"""


import os, sys, pickle, time
from itertools import combinations

from data_society import Society
from input_functions import get_project_dict, _Minimal_Sample_Number
#_Minimal_Sample_Number: least number of shared subjects is required to perform network association.

from pls2_network import pairNetwork



class HiCoNet:
    """
    Hirarchical community network.
    Need a project dictionary to initiate, either from local project.yaml file or from web input.

    Each data type is a society, where local communities are identified.
    HiCoNet is based on PLS association cross data types.
    Time points are not explicitly used but it's users' decision to place into project dictionary.

    Return network, and supporting database, traceback for data presentation.
    
    
    The pairNetworks can be combined.

    

    """
    def __init__(self, dict_project_definition):
        
        self.dict_project_definition = dict_project_definition
        
        # get self.societies and self.society_dict
        self.societies = []
        self.get_Societies()
        
        

    def run_hiconet(self):
        # verify association_dict and update if needed
        self.get_association_instructions()
        
        # heavy lifting - PLS2 regression and permutation for each association
        self.networks = []
        self.run_community_associations()
        
        self.make_outdir()
        
        self.write_networks()
        self.write_top_networks()
        self.write_societies()
        self.pickle_hiconet()
        
    def get_association_instructions(self):
        """ 
        to get associations to compute: check user's spec or infer ab initio.
        if no valid instruction from user input, do ab initio inference.
        
        and ?? Complicated scheme
        
        allow cross-timepoint too...
        
        Update self.dict_project_definition['associations']
        """
        # if not self.dict_project_definition.has_key('associations')
        if not 'associations' in self.dict_project_definition:
            self.dict_project_definition['associations'] = []
            for sc1,sc2 in combinations( self.societies, 2 ):
                print("Inferring association instruction", sc1.name, sc2.name)
                self.dict_project_definition['associations'] += self.infer_association_dicts(sc1, sc2)

        else:
            print("Checking association instruction...")
            self.dict_project_definition['associations'] = [self.check_update_assoc_dict(A, 
                                                            self.society_dict[A['society1']], self.society_dict[A['society2']]) 
                                                            for A in self.dict_project_definition['associations']]
            
    def check_update_assoc_dict(self, A, sc1, sc2):
        """
        check/complete user's spec
        """
        # if user specified both lists of observation_list_society, use those
        if A['observation_list_society1'] and A['observation_list_society2']:
            # check equal length
            if len(A['observation_list_society1']) != len(A['observation_list_society2']):
                # define custom error later
                raise ValueError('Need equal numbers of observations to compute association.')
            else:
                # check same subject. Issue warning if mismatched. Don't force by overwriting user's intend
                if sc1.get_subjects_by_observation_list() != sc2.get_subjects_by_observation_list():
                    print("Warning - mismatched subjects in user specified observation_lists")
                
        # if IDs are specified (subjects are not necessarily corresponding to observation IDs) for one list, find out the other
        # this does not limit to single time point
        elif A['observation_list_society1']:
            # find matched observation_list_society2
            subjects = sc1.get_subjects_by_observation_list()
            timepoints = sc1.get_timepoints_by_observation_list()
            A['observation_list_society2'] = sc2.get_observation_list(subjects, timepoints)
                
        elif A['observation_list_society2']:        
            # find matched observation_list_society1
            subjects = sc2.get_subjects_by_observation_list()
            timepoints = sc2.get_timepoints_by_observation_list()
            A['observation_list_society2'] = sc1.get_observation_list(subjects, timepoints)
            
        else:
            # If neither observation_lists is specified
            A = {}
            
        return A
        
        
    def infer_association_dicts(self, sc1, sc2):
        """
        infer from ObservationAnnotation in both societies, crossing all time points.
        
        associations:
        - name: gene_ab
          society1: genes
          observation_list_society1: []
          society2: antibody
          observation_list_society2: []
        
        Returns
        -------
        list_associations, [{}, ...]
        
        """
        
        def make_dict_(L1, L2, sc1, sc2, tp1, tp2):
            # based on validated lists [(subj, observation_ID), ...]; matching subjects
            A = {}
            dict2 = dict(L2)
            L1_use = [x for x in L1 if x[0] in dict2.keys()]
            A['observation_list_society1'] = [x[1] for x in L1_use]
            A['subjects'] = [x[0] for x in L1_use]
            A['observation_list_society2'] = [dict2[subj] for subj in A['subjects']]
            A['society1'] = sc1.name
            A['society2'] = sc2.name
            A['name'] = sc1.name + '_' + sc2.name
            A['timepoint1'] = tp1
            A['timepoint2'] = tp2
            return A
        
        list_associations = []
        for T1 in sc1.timepoints:
            for T2 in sc2.timepoints:
                # get obser in T
                L1 = [(sc1.dict_obser_subj[obser], obser) 
                      for obser in sc1.DataMatrix.columns if sc1.dict_obser_time[obser]==T1]
                L2 = [(sc2.dict_obser_subj[obser], obser) 
                      for obser in sc2.DataMatrix.columns if sc2.dict_obser_time[obser]==T2]                
                list_associations.append( make_dict_(L1, L2, sc1, sc2, T1, T2) )
            
        return list_associations
            

    def run_community_associations(self):
        """
        If user defined lists of observartions to use, use them.
        Otherwise, infer from sample lists.
        
        association_dict contains how the association is to be computed: matching samples etc.
        
        pn.network_edges = [( g, m, PLSscore, p-value ), ...]
        """
        for A in self.dict_project_definition['associations']:
            if len(A['subjects']) > _Minimal_Sample_Number:
                #print(A)
                pn = pairNetwork(A, self.society_dict[A['society1']], self.society_dict[A['society2']])
                self.networks.append( pn )


    def get_Societies(self):
        """
        Generate a list and dictionary of societies.
        """
        for sd in self.dict_project_definition['societies']:
            if sd['file_data_matrix']:
                print("Working on {}".format(sd))
                self.societies.append(Society(sd, self.dict_project_definition['workdir']))

        self.society_dict = {}
        for sc in self.societies:
            #
            # to add verification that society name is unique
            #
            self.society_dict[sc.name] = sc

    def parse_yaml_design(proj_dict):
        """
        Deduce the instruction on project,
        what data types, how to integrate.
        If not specified in valid yaml, 
        do default, to get PLS association btw all datatype/timepoints, 
        including diff to baseline (proper time course analysis for future development).
        
        """
        pass

    def make_outdir(self):
        time_stamp = str(time.time())
        self.outdir = os.path.join(self.dict_project_definition['workdir'], 
                                   'output_hiconet_'+time_stamp)
        os.mkdir(self.outdir)
        
    def write_networks(self):
        """
        
        
        
        
        create local DIR, mummichog style
        
        
        
        
        
        pn.network_edges = [( g, m, PLSscore, p-value ), ...]
        """
        
        for pn in self.networks:
            s = "datatype_1\tdatatype_2\tcommunity1_number\tcommunity2_number\tPLS_score\tp-value\n"
            for e in pn.network_edges:
                s += '\t'.join( [pn.dict['society1'], pn.dict['society2']] + [str(x) for x in e]) + '\n'

            with open(os.path.join(self.outdir, pn.name + "_network.txt"), "w") as O:
                O.write(s)


    def write_top_networks(self):
        all = []
        for pn in self.networks:
            for e in pn.network_edges:
                all += [[pn.name] + list(e)]
        
        def sort2(val): return val[3]
        all.sort(key = sort2, reverse = True)
        top = [x for x in all if x[4] < 0.05 and x[3] > 0]
        if len(top) < 10:
            top = all[:10]
        
        #print(top)
        
        with open(os.path.join(self.outdir, pn.name + "_network.txt"), "w") as O:
            O.write( '\n'.join(['\t'.join([str(ii) for ii in x]) for x in top]) )


    def write_societies(self):
        """
        Writes out each society, and its communities
        """
        for sc in self.societies:
            with open(os.path.join(self.outdir, sc.name + "_communities.txt"), "w") as O:
                O.write(str(sc.Communities))

    def pickle_hiconet(self):
        """
        write Python pickle object for entire HiCoNet instance
        

        for k,v in self.comm_member_dict.items():
            # cluster number : row numbers
            self.dataframe.iloc[v, :].to_csv( OUTDIR + self.datatype + "_clus_%d.txt" %k, sep="\t")
        
        """
        pass


class DeltaHiCoNet(HiCoNet):
    """
    This copies HiCoNet, but uses differences btw time points for PLS regression.
    
    to do.
    
    """
    def make_delta_societies(self, old_hiconet):
        new_societies = []
        for sc in self.societies:
            sc.time_points = []
            sc.DataMatrix = []
            
        self.societies = new_societies





if __name__ == '__main__':
    '''
    Run analysis

    DH = DeltaHiCoNet()
    DH.make_delta_societies(H)
    DH.run_hiconet()


    '''

    dir = sys.argv[1]
    proj_dict = get_project_dict(_source='local', _dir=dir)     #_Data_Directory)
    H = HiCoNet(proj_dict)
    H.run_hiconet()
    

    
    
