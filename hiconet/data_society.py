
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

from input_functions import read_input_tables, fuzzy_index_2L, data_wrangler, \
                            common_observation_IDs, common_subject_IDs, common_timepoint_labels, common_treatment_labels
from community_detection import hierachical_clustering

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

        """
        [
        self.raw_DataMatrix,
        self.ObservationAnnotation,
        self.FeatureAnnotation,
        self.unstructured,
        ] = read_input_tables(_dir, dict_society)

        # dict feature annotation
        self.get_feature_annotation()
        # subj _ time point _ treatment group
        self.get_observation_annotation()

        # clean up DataMatrix
        self.cleanup()
        
        # delta should be considered as a new society
        #self.generate_delta()

        # match annotation
        # sample_map_check
        # summary, issues, stats
        


    def get_feature_annotation(self):
        # ??
        #return dict(self.FeatureAnnotation)
        pass

    def get_observation_annotation(self):
        """
        To get a data structure 
        self.ObservationAnnotation_list = [[observation_ID, subj_ID, time point, treatment], ...] 
        
        This will be used to match btw socieities in constructing the association dicts
        
        `treatment` is unclear in ImmPort ontology?
        Not sure how to deal with missing index in pandas dataframe; manual index instead. 
        Use '' if not found.


        ObservationAnnotation may use different terms. Thus includign a list of common synonyms and variants.

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
        Filtering (both observations and features) and imputation
        
        Pass through for now - yet to implement data_wrangler.
        
        """
        df, message = data_wrangler(self.raw_DataMatrix)
        self.track_messages.append(message)
        if message:
            self.DataMatrix = df
        else:
            self.DataMatrix = self.raw_DataMatrix




    def get_communities(self, method='hcl'):
        """
        Testing stage, using hcl only
        
        option to designate minimal cluster number/size?

        self.Communities is a dictionary {community_ID: [feature_ID, ...], ...}
        """
        available_methods = ['hcl',
                            'lcms_hcl',
                            'spec-Girvan-Newman',
                            'Louvain',
                            # 'leidenalg',
                            'test',
                            ]
        if method not in available_methods:
            raise ValueError('Provide a valid method, one of {}.'.format(available_methods))
        if method == 'hcl':
            # This return format will change -
            self.Clus, self.Communities = hierachical_clustering(self.DataMatrix)
        if method == 'Louvain':
            pass
        elif method == '':
            pass


    def get_comm_member_dict(self):
        """
        This is already in hierachical_clustering
        
        Not used now, but should be moved out of hierachical_clustering
        
        return {community_ID: [feature_ID, ...], ...}
        """
        d = {}
        for k,v in self.Communities.items():      # {feature_ID: community_ID}
            if v in d:
                d[v].append(k)
            else:
                d[v] = [k]
        return d



    def write_communities(self):

        # ?? not to use for now
        for k,v in self.comm_member_dict.items():
            # cluster number : row numbers
            self.DataMatrix.iloc[v, :].to_csv( OUTDIR + self.datatype + "_clus_%d.txt" %k, sep="\t")



    def __generate_delta(self):
        """
        This expands the DataMatrix and ObservationAnnotation to add time differences btw observations.
        E.g. antibody increase from baseline.

        # delta should be considered as a new society

        """
        def make_observation_id(id1, id2): return '_minus_'.join((str(id2), str(id1)))

        pass



