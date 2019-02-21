"""Hieracrchical Community Network - data I/O

from project.yaml or web input (JSON) to derive a list of dictionaries, 
one for each data type as dict_files_designation

---- YAML example ----
societies:
- name: genes
  datatype: transcriptomics
  file_data_matrix: ''
  file_feature_annotation: ''
  file_observation_annotation: ''
  file_unstructured: ''

- name: cells
  datatype: facs
  file_data_matrix: ''
  file_feature_annotation: ''
  file_observation_annotation: ''
  file_unstructured: ''

# use load, not load_all
>>> j = yaml.load_all(open('project.yaml').read())
>>> for i in j:
...   print(i)
... 
{'project': 'HiCoNet_test_80', 'source_id': 'SDY80', 'retrieved_data': '2018-0-0', 'data_from': 'ImmPort', 
'societies': [{'file_feature_annotation': '', 'file_observation_annotation': '', 'file_unstructured': '', 'name': 'genes', 'datatype': 'transcriptomics', 'file_data_matrix': 'subjLabel_BTM_SDY80_pH1N1_2009_Transcriptomics.txt'}, 
{'file_feature_annotation': '', 'file_observation_annotation': '', 'file_unstructured': '', 'name': 'cells', 'datatype': 'facs', 'file_data_matrix': ''}, 
{'file_feature_annotation': '', 'file_observation_annotation': '', 'file_unstructured': '', 'name': 'antibody', 'datatype': 'HAI', 'file_data_matrix': ''}, 
{'file_feature_annotation': '', 'file_observation_annotation': '', 'file_unstructured': '', 'name': 'metabolites', 'datatype': 'metabolomics', 'file_data_matrix': ''}]}
>>> 


workdir is also passed to proj_dict after reading yaml.
"""
import os, yaml
import pandas as pd
from fuzzywuzzy import process as fuzzyfind

# Will need better management of parameters
#_Data_Directory = './datasets/SDY80/'
_Minimal_Sample_Number = 10
NUM_PERMUTATION = 20    #200

common_observation_IDs = ['BIOSAMPLE_ACCESSION', 'observation_ID', 'sample']
common_subject_IDs = ['SUBJECT_ACCESSION', 'Participant ID', 'SUBJECT_ID', 'subject', 'Participant', 'patient', 'patient_ID', 'volunteer']
common_timepoint_labels = ['STUDY_TIME_COLLECTED', 'Study Time Collected', 'time point']
common_treatment_labels = ['treatment', 'infected']


def read_input_tables(_dir, dict_society):
    """Process input files into dataframes/dict for HiCoNet

    Parameters
    ----------
    dict_society:
        {'name': 'genes', 
        'datatype': 'transcriptomics', 
        'file_data_matrix': 'subjLabel_BTM_SDY80_pH1N1_2009_Transcriptomics.txt',
        'file_feature_annotation': '', 
        'file_observation_annotation': '', 
        'file_unstructured': '', }
    All tables have to comply with format: data matrix is M x N, features as col 0.
    Annotation tables use col 0 as common ID
    unstructured is optional.

    Returns
    -------
    unprocessed, unlinked tables:
    [DataMatrix, ObservationAnnotation, FeatureAnnotation, unstructured]
    """
    DataMatrix = pd.read_csv(
        os.path.join(_dir, dict_society['file_data_matrix']), sep='\t', index_col=0)  # index_col=0 specifies format for the feature column.
                                                        # pandas can use 1st col as feature index if the 1st row misses 1st cell, 
                                                        # but the R convention should be avoided here for more general applications
    ObservationAnnotation = pd.read_csv(
        os.path.join(_dir, dict_society['file_observation_annotation']), sep='\t')    # not using index_col=0

    FeatureAnnotation = {}
    if dict_society['file_feature_annotation']:
        FeatureAnnotation = pd.read_csv(os.path.join(_dir, dict_society['file_feature_annotation']), sep='\t')

    unstructured = {}
    if dict_society['file_unstructured']:
        unstructured = {}

    return [DataMatrix, ObservationAnnotation, FeatureAnnotation, unstructured]



def get_project_dict(_source='local', _dir='.'):
    """If local, look for file 'project.yaml'; (else, get from web?)

    Returns
    -------
    project dictionary
    """
    if _source == 'local':
        files = os.listdir(_dir)
        if 'project.yaml' in files:
            d = read_yaml_file(os.path.join(_dir, 'project.yaml'))
            d['workdir'] = _dir
            return d
        else:
            raise Exception("'project.yaml' not found.")

def read_yaml_file(yaml_file):
    return yaml.load(open(yaml_file).read())



def data_wrangler(df):
    """
    Perform data cleaning and imputation.
    Filter by missing values, both samples (observations) and features.
    Impute the rest of missing values.
    track steps and warnings. 
    
    to do.

    Returns
    -------
    (cleaned up dataframe, process record)

    
    """
    return (df, '')


def fuzzy_index(target, L):
    """
    Find the index of target in list L.
    Allowing fuzzy match, accommodating typos and cases etc.
    
    Returns
    -------
    [index, top matched target in L]
    """
    t = fuzzyfind.extractOne(target, L)
    # limited to good matches. t = (matched, score)
    if t[1] > 80:
        return [L.index(t[0]), t[0]]
    else:
        return [None, None]

def fuzzy_index_2L(targetList, L):
    """
    Find the index of target in list L. targetList is used to include common variants and synonyms.
    Allowing fuzzy match, accommodating typos and cases etc.
    
    Returns
    -------
    [index, top matched target in L]
    """
    t = (None, 0)
    for target in targetList:
        t1 = fuzzyfind.extractOne(target, L)
        if t1[1] > t[1]:
            t = t1
            
    if t[1] > 80:
        return [L.index(t[0]), t[0]]
    else:
        return [None, None]


def data_indexer(df):
    """Data indexing by subject and timepoint.
    Will remove all quotes.
    
    Not used now.

    In [34]: bigdict = {}
    
    In [35]: for line in w[1:]:
        ...:     # construct index subj_timepoint_cellpopulation
        ...:     a = line.split('\t')
        ...:     k = a[1].replace('"', '').replace('.372', '') + '_' + a[6] + '_' + a[11].replace('"', '')
        ...:     bigdict[k] = a[8]
        ...:     
    
    In [36]: bigdict.items()[:4]
    Out[36]: 
    [('SUB135758_90_CD8+ T cell', '"510.6118469238281"'),
     ('SUB135773_28_Tfh17 CD4+ T cell', '"21.569305419921875"'),
     ('SUB135772_28_CD40- BDCA3+ myeloid dendritic cell', '"0.37419480085372925"'),
     ('SUB135772_0_HLA-DR+ Th1 CD4+ T cell', '"0.06463897228240967"')]
    
    In [37]: # construct wanted list
    
    In [38]: s = 'CellPopulation\t' + '\t'.join(populations) + '\n'
    
    In [39]: row = [x+"_0_" for x in subjects] + [x+"_1_" for x in subjects] + [x+"_7_" for x in subjects] + [x+"_28_" for
        ...:  x in subjects] + [x+"_90_" for x in subjects]
    
    In [40]: len(row)
    Out[40]: 85
    
    In [41]: row[:4]
    Out[41]: ['SUB135758_0_', 'SUB135759_0_', 'SUB135760_0_', 'SUB135761_0_']
    
    In [42]: s = 'CellPopulation\t' + '\t'.join(row) + '\n'
    
    In [43]: for p in populations:
        ...:     s += p + '\t' + '\t'.join([bigdict.get(x+p, 'NA') for x in row]) + '\n'
        ...:     
    
    In [44]: with open('formatted_SDY372_fcs_analyzed_result.txt', 'w') as O:
        ...:     O.write(s)
        ...:     
    
    In [45]: # quotes included in the output above
    


    """
    value_dict = {}
    return value_dict



