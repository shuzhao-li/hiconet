
"""
Hieracrchical Community Network
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

"""
import pandas as pd
try:
    from pandas.compat import StringIO
except ImportError:
    from io import StringIO

from .input_functions import read_input_tables, fuzzy_index_2L, data_wrangler, \
                            common_observation_IDs, common_subject_IDs, common_timepoint_labels, common_treatment_labels, \
                            auto_BTM_conversion, gene_2_btm
from .community_detection import hierachical_clustering, leiden_find_communities, hierachical_clustering_lcms

class Society:
    """
    Each data type = a Society,
    From which local feature communities are identified.
    The communities are defined by a dictionary, {feature_ID: community_ID}.


    To define time points, sample-subject matching, we parse the dataframe of ObservationAnnotation
    into a list of [[observation_ID, subj, time point, treatment], ...]

    Time points or treatments are used in slicing data for further analysis - Think those as sample groups.

    """
    def __init__(self, dict_society, _dir):
        """
        Initiation from files.
        Given the size of these files, the web version will aslo take files and store first.
        Not a good idea to read file directly from browser. Will add compatibility of compressed formats.

        dict_society:
        {'name': 'genes',
        'datatype': 'transcriptomics',
        'file_data_matrix': 'subjLabel_BTM_SDY80_pH1N1_2009_Transcriptomics.txt',
        'file_feature_annotation': '',
        'file_observation_annotation': '',
        'file_unstructured': '', }
        
        self.Communities is a dictionary {community_ID: [feature_ID, ...], ...}
        self.Clus is [feature_community_number, ...]

        """
        self.name = dict_society['name']
        self.datatype = dict_society['datatype']
        self.graph = []
        self.Communities = {}
        self.track_messages = []
        
        # these are overall subjects and time points, not necessarily in data matrix
        self.subjects = []
        self.timepoints = []
        
        self.populate_data_tables(dict_society, _dir)
        self.get_communities()


    def populate_data_tables(self, dict_society, _dir):
        """
        Get input tables, also match annotation, check stats
        #[DataMatrix, ObservationAnnotation, FeatureAnnotation, unstructured,] = read_input_tables(_Data_Directory, dict_society)

        # summary, issues, stats
        """
        [
        self.raw_DataMatrix,
        self.ObservationAnnotation,
        self.FeatureAnnotation,
        self.unstructured,
        ] = read_input_tables(_dir, dict_society)

        # clean up DataMatrix
        self.cleanup()
        if auto_BTM_conversion:
            self.auto_BTM_convert()
        
        # dict feature annotation
        self.get_feature_annotation()
        # subj _ time point _ treatment group
        self.get_observation_annotation()

        
    def auto_BTM_convert(self):
        """Convert gene table to BTM (Blood Transcription Modules) activity table.
        
        ?? need update project dictionary?
        """
        # check is_geneTable; only convert if > 3000 features
        
        if self.datatype == 'transcriptomics' and self.DataMatrix.shape[0] > 3000:
            print("Converting transcriptomics data to BTM modules..")
            self.DataMatrix = gene_2_btm(self.DataMatrix)
        

    def get_feature_annotation(self):
        # ??
        # temporary hack for feature IDs; will do more proper annotation
        #
        self.feature_member_annotation = list(self.DataMatrix.index)

    def get_observation_annotation(self):
        """
        To extract a new data structure from self.ObservationAnnotation
        self.ObservationAnnotation_list = [[observation_ID, subj_ID, time point, treatment], ...] 
        
        This will be used to match btw socieities in constructing the association dicts
        
        `treatment` is unclear in ImmPort ontology?
        Not sure how to deal with missing index in pandas dataframe; manual index instead. 
        Use '' if not found.


        ObservationAnnotation may use different terms. Thus including a list of common synonyms and variants.

        """
        self.ObservationAnnotation_list = []
        # potential to improve the lists of common synonyms and variants
        wanted = [fuzzy_index_2L(T, list(self.ObservationAnnotation.columns))[0] 
                        for T in [common_observation_IDs, common_subject_IDs, common_timepoint_labels, common_treatment_labels]]
        
        # likely missing 'treatment' or something
        print(wanted)
        
        for L in self.ObservationAnnotation.values:
            tmp = [''] * len(wanted)
            for ii in range(len(wanted)):
                if isinstance(wanted[ii], int):
                    tmp[ii] = L[wanted[ii]]
            self.ObservationAnnotation_list.append(tmp)

        #print("self.ObservationAnnotation_list", self.ObservationAnnotation_list)

        self.subjects = set([ L[1] for L in self.ObservationAnnotation_list ])
        self.timepoints = set([ L[2] for L in self.ObservationAnnotation_list ])
        
        print("self.timepoints", self.timepoints)
        
        self.make_dict_obser_subj()
        self.make_dict_obser_time()
        

    def make_dict_obser_subj(self):
        """
        Retrieve subjects based on 
        self.ObservationAnnotation_list = [[observation_ID, subj_ID, time point, treatment], ...]
        """
        d = {}
        for L in self.ObservationAnnotation_list: d[L[0]] = L[1]
        self.dict_obser_subj = d

    def get_subjects_by_observation_list(self, L):
        return [self.dict_obser_subj.get(x, None) for x in L]
        
    def make_dict_obser_time(self):
        """
        Retrieve time points based on 
        self.ObservationAnnotation_list = [[observation_ID, subj_ID, time point, treatment], ...]
        """
        d = {}
        for L in self.ObservationAnnotation_list: d[L[0]] = L[2]
        self.dict_obser_time = d
        
        
    def get_timepoints_by_observation_list(self, L):
        return [self.dict_obser_time.get(x, None) for x in L]
    
    def get_observation_list(self, subjects, timepoints):
        """
        Retrieve observations based on 
        self.ObservationAnnotation_list = [[observation_ID, subj_ID, time point, treatment], ...]
        """
        d = {}
        for L in self.ObservationAnnotation_list: d[(L[1], L[2])] = L[0]
        return [d.get(x, None) for x in zip(subjects, timepoints)]



    def cleanup(self):
        """
        Not really cleaning now -
        
        Filtering (both observations and features) and imputation
        
        Pass through for now - yet to implement data_wrangler.
        
        """
        df, message = data_wrangler(self.raw_DataMatrix)
        self.track_messages.append(message)
        self.DataMatrix = df
        if message:
            #print(df.values)
            print(message)


    def get_communities(self, method=''):
        """
        Default is leiden method for commuity detection.
        If not specified, lcms_hcl is default for LC-MS metabolomics.
        
        option to designate minimal cluster number/size, and optimal number of clusters.

        self.Communities is a dictionary {community_ID: [feature_ID, ...], ...}
        
        TO-DO add community annotation function
        
        """
        available_methods = ['hcl',
                            'lcms_hcl',
                            #'spec-Girvan-Newman',
                            #'Louvain',
                            'leiden',
                            #'test',
                            ]
        if not method:
            if self.datatype == 'metabolomics': method = 'lcms_hcl'     # this needs further specifications for more data formats
            else: method = 'leiden'
        
        # dispatch 
        if method == 'leiden':
            self.Clus, self.Communities = leiden_find_communities(self.DataMatrix)
        elif method == 'lcms_hcl':
            # also using information, e.g. retention time, in FeatureAnnotation
            print("Using lcms_hcl for LC-MS metabolomics.")
            self.Clus, self.Communities = hierachical_clustering_lcms(self.DataMatrix, self.FeatureAnnotation)
        elif method == 'hcl':
            # is this a good return format?
            self.Clus, self.Communities = hierachical_clustering(self.DataMatrix)
        elif method not in available_methods:
            raise ValueError('Provide a valid method, one of {}.'.format(available_methods))

        #self.annotate_communities()


    def export_communities_table(self):
        """
        TO-DO: improve community annotation
        """
        s = "community_number\tfeature_number\tfeature_name\n"
        for k,v in self.Communities.items():
            for member in v:
                s += '\t'.join([str(x) for x in 
                                [k, member, self.feature_member_annotation[member]]]) + '\n'
        return s



class webSociety(Society):
    """Simplified for web version.
    self.get_communities() is run upon init.
    """
    
    def populate_data_tables(self, dict_society, _dir=''):
        """
        Get input tables, also match annotation, check stats
        #[DataMatrix, ObservationAnnotation, FeatureAnnotation, unstructured,] = read_input_tables(_Data_Directory, dict_society)
        dict_society:
            {'name': 'genes', 
            'datatype': 'transcriptomics', 
            'file_data_matrix': input_str,
            }
        All tables have to comply with format: data matrix is M x N, features as col 0.
        Annotation tables use col 0 as common ID
        
        For web input, pre-matched data columns are required, thus simplified data matching.
        """
        self.raw_DataMatrix = pd.read_csv(
            StringIO(dict_society['file_data_matrix']), sep='\t', index_col=0)
        # pre-matched data columns are required for web I/O
        # ObservationAnnotation has to be constructed in a way to be compatible with default data structure 
        self.ObservationAnnotation = {}
        # FeatureAnnotation may be necessary, e.g. LC-MS data are grouped using rtime in FeatureAnnotation
        
        self.FeatureAnnotation = {}
        if dict_society['file_feature_annotation']:
            self.FeatureAnnotation = pd.read_csv(
                StringIO(dict_society['file_feature_annotation']), sep='\t', index_col=0)

        
        # clean up DataMatrix
        self.cleanup()
        if auto_BTM_conversion:
            self.auto_BTM_convert()
                
        # dict feature annotation
        self.get_feature_annotation()
        # subj _ time point _ treatment group
        self.get_observation_annotation()
        
    def get_observation_annotation(self):
        """
        simplified for web I/O, which requires prematched samples.
        
        self.ObservationAnnotation_list = [[observation_ID, subj_ID, time point, treatment], ...] 
        
        This will be used to match btw socieities in constructing the association dicts
        
        Set for now: self.timepoints = [0]
        
        ------
        dictionary {'user_supplied_0': 'col0', 'user_supplied_1': 'col1', ...}
        """
        self.ObservationAnnotation_list = []
        ii = 0
        for observation_ID in self.raw_DataMatrix.columns:
            ii += 1
            self.ObservationAnnotation_list.append( [observation_ID, 'col'+str(ii), 0, 'x'] )

        #print("self.ObservationAnnotation_list", self.ObservationAnnotation_list)
        self.subjects = set([ L[1] for L in self.ObservationAnnotation_list ])
        self.timepoints = set([ L[2] for L in self.ObservationAnnotation_list ])
        
        #print("self.timepoints", self.timepoints)
        self.make_dict_obser_subj()
        self.make_dict_obser_time()






