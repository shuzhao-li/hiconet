# Hieracrchical Community Network
for integration of multiple data types collected from a common group of subjects.

## Terminology:
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


    Graph: a graph/network for relationships in the data (e.g. used in loom format, loompy.org).
    Community: a group of features within a society that share a similar pattern.



    Reusing
    -------
    The data structure is aligned with anndata (https://github.com/theislab/anndata) but transposed.
    For future consideration, e.g. in loom format, one loom file = a Society; as meta data can differ for different data types.
    Same goes for anndata.


## requires
    'numpy',
    'scipy',
    'pandas',
    'sklearn',
    'leidenalg',
    'igraph',
    'fuzzywuzzy',

    Note: python-igraph requires the C library igraph  
    https://stackoverflow.com/questions/45667147/install-python-igraph-on-mac


## Test run
➜  hiconet python hiconet/HiCoNet.py hiconet/datasets/SDY80

[Output like -]
/usr/local/lib/python2.7/site-packages/fuzzywuzzy/fuzz.py:11: UserWarning: Using slow pure-python SequenceMatcher. Install python-Levenshtein to remove this warning
  warnings.warn('Using slow pure-python SequenceMatcher. Install python-Levenshtein to remove this warning')
Working on {'file_feature_annotation': '', 'file_observation_annotation': 'biosample.txt', 'file_unstructured': '', 'name': 'genes', 'datatype': 'transcriptomics', 'file_data_matrix': 'BTMactivity_SDY80.txt'}
[0, 9, 5, None]
('self.timepoints', set([0, 1, 70, 7, -7, 28]))
((346, 292), (59685,))
[12  7  7  7 10  7 10  6  6  6 10  6  6  6  6  6  6  6  6  6 10  6 10 10
  1  6  2  1  2  2  1  6 10  6  6 10 10 10  1 13  1  6 10  7  7  7  7  1
  1 12  5  4  4 10 12 10  7  7  7 13 13 10 12 11 12  8  8  8  8  8  8  8
  8  8 10 11 12  2  8 10 10 10  4 13  8 10  7  7  7  7  7 10 12 12  1  8
  6  5  5  5  5  5 10  6 10 11  8 10 10 10 10  5 10 10  5  2  2  2  9  9
  8 10  1 10 13 13  5  7  7 10  3  3  3 10  9 13  4  7 12 10  9 11  7  5
  8 11 12 13 10  7  7  7 12 12  7  2 10  2 12 10 10  7  5  7  7  7  8  3
  5  6 11 10  4  4  7  7 10  7 10 13 10 11 10 10 10 13  7  1 10 10 10  8
  6  7 12  9  7  1 13  3  3  2 11 10  7  7  7  7  7  8  3  3 10 10 12  7
  4  3 10 11 10  3  3  7 13  7  9  9  9  3  7  7 10  5  6  2  9  9 11 12
  3  7  7 10  9 13  7  9 10  4  7  9 12  7  3  7  7  4  4  9  3  3  6  4
  7  7  7  3  7  7 10  7 12  2  8 11  4  7 11  7 10 11 10  9  7  3  4  4
  3  7  7  2 10  7  3  6  3  3  8  6  7  4  6  7  9  9  2  7  9  4  4  7
  3  3  6  4  7  8  4  4  4  6  7  2  2  3  8  7  6  3  2  7  7  6  1  2
  9  7 10 10  2  1  9  9 10 10]
('number of features: ', 346)
('number of clusters: ', 13)
[(1, [24, 27, 30, 38, 40, 47, 48, 94, 122, 187, 197, 334, 341]), (2, [26, 28, 29, 77, 115, 116, 117, 155, 157, 201, 235, 273, 291, 306, 323, 324, 330, 335, 340]), (3, [130, 131, 132, 167, 199, 200, 210, 211, 217, 221, 222, 229, 240, 254, 260, 261, 267, 285, 288, 294, 296, 297, 312, 313, 325, 329])]
Working on {'file_feature_annotation': '', 'file_observation_annotation': 'neut_ab_titer__observation_annotation.txt', 'file_unstructured': '', 'name': 'antibody', 'datatype': 'antibody', 'file_data_matrix': 'neut_ab_titer__data_matrix.txt'}
[0, 1, 6, None]
('self.timepoints', set([0.0, -7.0, 28.0, 70.0, 7.0]))
((4, 240), (6,))
[1 1 1 1]
('number of features: ', 4)
('number of clusters: ', 1)






