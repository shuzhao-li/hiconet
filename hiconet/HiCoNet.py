# Copyright 2018-2019 Shuzhao Li. All Rights Reserved.
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

See README for details.

"""


__version__ = "0.5.3"
__updated__ = "2023-02-13"


import os
import sys
import pickle
import time
import json
from itertools import combinations

from .data_society import Society, webSociety
from .input_functions import get_project_dict, _Minimal_Sample_Number

from .pls2_network import pairNetwork
from .html_visual import make_js_from_network, make_html_page


class HiCoNet:
    """
    Hirarchical community network.
    Need a project dictionary to initiate, either from local project.yaml file or from web input.

    Each data type is a society, where local communities are identified.
    HiCoNet is based on PLS association cross data types. The pairNetworks can be combined.
    
    For web server, we can start a simplified input format without considering time points.
    
    Return
    ------
    One community network, 
    List of communities,
    Definition of each community.
    
    Output formats can be
    JSON, pickle and local file writes.

    """
    def __init__(self, dict_project_definition):
        
        self.dict_project_definition = dict_project_definition
        
        # get self.societies and self.society_dict
        self.societies = []
        self.get_Societies()
        self.networks = []
        self.network_dict = {}
        

    def run_hiconet(self):
        """
        Local run routine
        """
        # verify association_dict and update if needed
        self.get_association_instructions()
        
        # heavy lifting - PLS2 regression and permutation for each association
        self.run_community_associations()
        
        self.make_outdir()
        
        self.sort_combined_network()
        self.write_networks()
        self.write_societies()
        self.pickle_hiconet()
        
    def get_association_instructions(self):
        """ 
        to get associations to compute: check user's spec or infer ab initio.
        if no valid instruction from user input, do ab initio inference.
        
        """
        
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
            # L1_use defines common subjects
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
                # get obser in the specified time point
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
                print("\n############################\n\n")
                print(A)
                pn = pairNetwork(A, self.society_dict[A['society1']], self.society_dict[A['society2']])
                self.networks.append( pn )
                self.network_dict[pn.name] = pn
                 
            else:
                print("Dropped %s, fewer samples than minimal requirement." %A['name'])


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

    def parse_yaml_design(self, proj_dict):
        """
        Deduce the instruction on project,
        what data types, how to integrate.
        If not specified in valid yaml, 
        do default, to get PLS association btw all datatype/timepoints, 
        including diff to baseline (proper time course analysis for future development).
        
        """
        pass

    def make_outdir(self):
        """
        create local DIR, mummichog style
        """
        time_stamp = str(time.time())
        self.outdir = os.path.join(self.dict_project_definition['workdir'], 
                                   'output_hiconet_'+time_stamp)
        self.outdir_individual_networks = os.path.join(self.outdir, 'individual_networks')
        os.mkdir(self.outdir)
        os.mkdir(self.outdir_individual_networks)
        
    def write_networks(self):
        """
        Write all network edges, and a separate file for top networks.
        But the top networks use composite community IDs.
        """
        header = "datatype_1\tdatatype_2\tcommunity1_number\tcommunity2_number\tPLS_score\tp-value\n"
        for pn in self.networks:
            s = header
            for e in pn.network_edges:
                # pn.network_edges = [( g, m, PLSscore, p-value ), ...]
                s += '\t'.join( [pn.dict['society1'], pn.dict['society2']] + [str(x) for x in e]) + '\n'

            with open(os.path.join(self.outdir_individual_networks, pn.name + "_network.txt"), "w") as O:
                O.write(s)

        # self.combined_network is created in sort_combined_network
        top = [x for x in self.combined_network if x[4] < 0.05 and x[3] > 0]
        if len(top) < 10:
            top = self.combined_network[:10]
            
        header2 = "datatypes_timepoints\tcommunity1_number\tcommunity2_number\tPLS_score\tp-value\n"
        with open(os.path.join(self.outdir, "top_networks.txt"), "w") as O:
            O.write(header2 + '\n'.join(['\t'.join([str(ii) for ii in x]) for x in top]) )
            
            
    def sort_combined_network(self):
        """Sort combined networks by p value
        Use composite community_numbers (`sc.name`.`community_ID`) from this function, 
        so that later JSON export can have unique ID for each community
        """
        all = []
        for pn in self.networks:
            for e in pn.network_edges:
                # pn.network_edges = [( g, m, PLSscore, p-value ), ...]
                all += [[pn.name, '.'.join((pn.dict['society1'], str(e[0]), str(pn.dict['timepoint1']))),
                         '.'.join((pn.dict['society2'], str(e[1]), str(pn.dict['timepoint2']))), 
                         e[2], e[3]]]
        
        def sort2(val): return val[4]
        all.sort(key = sort2, reverse = False)
        self.combined_network = all

    def write_societies(self):
        """
        Writes out each society, and its communities
        """
        for sc in self.societies:
            with open(os.path.join(self.outdir, sc.name + "_communities.txt"), "w") as O:
                O.write(sc.export_communities_table())

    def pickle_hiconet(self):
        """
        write Python pickle object for entire HiCoNet instance.
        Pickle is unsafe and should not be used for data exchange.
        
        
        To-do:
        more controlled dump.
        
        """
        pickle.dump(self, open(os.path.join(self.outdir, "hiconet_data.pickle"), "wb"))
        print("HiCoNet dumped into hiconet_data.pickle.")
        

    def export_json(self, max_num_edges=10):
        """export
        a) network of communities, b) list of communities, c) community definition/members
        
        combined_network = [(name, type1_community.n_Ti, type2_community.m_Tj, PLSscore, p-value ), ...]
        Community_ID = `sc.name`.`community_ID`.`timepoint`
        
        The basic organizaiton is hierarchy of three levels:
        community network, communities, features. 
        One community network, List of communities, Definition and members of each community.
        
        To-do:
        Feature heatmaps generated from the society data should be specific to time point
        Community annotations.
        To list two tabs, nodes/edges on right panel.
        
        """
        selected_combined_network =  self.combined_network[ :max_num_edges ]
        nodes, elements = self.table_to_cy_js_elements(selected_combined_network)
        community_members, community_data = {}, {}

        for sc in self.societies:
            for c, L in sc.Communities.items():
                # Communities is a dictionary {community_ID: [feature_index, ...], ...}
                for T in sc.timepoints:
                    n = '.'.join((sc.name, str(c), str(T)))
                    # Only export data in use
                    if n in nodes:
                        M, D = [], []
                        nested_data_list = sc.DataMatrix.values[L, :]
                        for ii in range(len(L)):
                            feature_name = sc.feature_member_annotation[L[ii]]
                            M.append(feature_name)
                            D.append({'name': feature_name, 'data': list(nested_data_list[ii])})
                            
                        community_members[n] = M
                        # To make this specific to time point in future, using pn.dict['observation_list_society1']
                        community_data[n] = D
        
        nodes = list(nodes)
        nodes.sort()
        result = {
            # this is the network (nodes and edges), called elements in cytoscape.js
            'elements': elements,
            'list_communities': nodes,
            'community_members': community_members,
            'community_data': community_data,
            }
        
        return json.dumps(result)

    def table_to_cy_js_elements(self, selected_combined_network):
        """
        Convert the combined_network in table format to cytoscape.js JSON format:
        http://js.cytoscape.org/#notation/elements-json
        
        combined_network = [(name, community_type1, community_type2, PLSscore, p-value ), ...]
        already sorted by p-value.
        
        
        Issues:
        edge ID has to be unique
        community IDs need to be specific w/ time point, etc.
        
        
        
        Returns
        =======
        Set of nodes, list of elements
        """
        allnodes = []
        for x in selected_combined_network: allnodes += x[1: 3]
        allnodes = set(allnodes)
        nodes = []
        for n in allnodes:
            # node type (classes) is already included in community IDs
            nodes.append({'data': {'id': n}, 'classes': [n.split('.')[0]],
                })
        
        edges = []
        for x in selected_combined_network:
            edges.append({'data': {'id': x[1]+'_'+x[2],
                                   'source': x[1],
                                   'target': x[2],
                                   'weight': x[4],  # p-value
                                   },
                })
        # cytoscape.js elements as flat array of nodes and edges
        return allnodes, nodes + edges



class DeltaHiCoNet(HiCoNet):
    """
    This copies HiCoNet, but uses differences btw time points for PLS regression.
    
    Yet to implement.
    
    """
    def make_delta_societies(self, old_hiconet):
        new_societies = []
        for sc in self.societies:
            sc.time_points = []
            sc.DataMatrix = []
            
        self.societies = new_societies

    def __generate_delta(self):
        """
        Not used. Yet to implement.
        
        This expands the DataMatrix and ObservationAnnotation to add time differences btw observations.
        E.g. antibody increase from baseline.

        # delta should be considered as a new society

        """
        def make_observation_id(id1, id2): 
            return '_minus_'.join((str(id2), str(id1)))

        pass



class webHiCoNet(HiCoNet):
    """
    HiCoNet via web I/O.
    
    Allowing scaling down, e.g. simplifed design, taking two prematched input matrices and compute HiCoNet.
    i.e. get two society, do community detection, then compute PLS2 network.
    
    dict_project_definition will be given by web input.
    From the web form, for each data type, a data matrix file, and optional annotation files are uploaded.
    
    """
        
    def get_Societies(self):
        """
        Generate a list and dictionary of societies.
        
        Community detection is performed at initiation of webSociety.
        
        not needed: 'file_unstructured': ''
        """
        
        for sd in self.dict_project_definition['societies']:
            if sd['file_data_matrix']:
                self.societies.append(webSociety(sd, ''))

        self.society_dict = {}
        for sc in self.societies:
            self.society_dict[sc.name] = sc

    def run_hiconet(self):
        """
        web run routine
        """
        # heavy lifting - PLS2 regression and permutation for each association
        self.get_association_instructions()
        print(self.dict_project_definition['associations'])
        
        self.run_community_associations()
        self.sort_combined_network()
        
        # to update for web output
        #print(self.export_json())


# -------------------------------------------------------------------
#
# Separate routines for regular HiCoNet, DeltaHiConet and webHiCoNet.
#
#
def run_local_HiCoNet():
    #
    proj_dict = get_project_dict(_source='local', _dir=sys.argv[1])
    H = HiCoNet(proj_dict)
    H.run_hiconet()
    
    # export JSON to file
    out = open(os.path.join(H.outdir, "result.json"), 'w')
    out.write(H.export_json())
    out.close()
    
    # testing html and js, up to 50 edges. This will be updated by the full JSON
    make_html_page(os.path.join(H.outdir, "summary.html"),
                   make_js_from_network(H.combined_network[:50]))


def run_web_HiCoNet(json_dict):
    """
    Use web form input to update project dictionary
    e.g.
    name_Field1, name_Field2, dataType1, dataType2, data_loc1, data_loc2
    """
    # parameters to be obtained via web; below is emptyp template.
    proj_dict = {
        'project': 'web_project_', 
        #'source_id': 'web_input', 
        #'retrieved_data': '2019-0-0', 
        #'data_from': 'web_user', 
        'societies': [
        {'name': 'data1', 
         'datatype': '', 
         'file_data_matrix': '',
         'file_feature_annotation': '', 
         'file_observation_annotation': '', 
         }, 
        {'name': 'data2', 
         'datatype': '', 
         'file_data_matrix': '',
         'file_feature_annotation': '', 
         'file_observation_annotation': '', 
         }, 
                      ]}
    
    proj_dict.update(json_dict
        )
        
        
    H = webHiCoNet(proj_dict)
    H.run_hiconet()
    


# -------------------------------------------------------------------
#
# local main
#
#

if __name__ == '__main__':
    '''
    Run analysis
    Separate routines for regular HiCoNet, DeltaHiConet and webHiCoNet.
    '''

    # regular HiCoNet from local dir
    run_local_HiCoNet()
    
    
    #run_web_HiCoNet()
    
    
