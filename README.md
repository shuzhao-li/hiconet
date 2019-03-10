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
    'PyYAML'
    'numpy',
    'scipy',
    'pandas',
    'sklearn',
    'leidenalg',
    'scanpy',
    'igraph',
    'fuzzywuzzy',

    Note: python-igraph requires the C library igraph  
    https://stackoverflow.com/questions/45667147/install-python-igraph-on-mac
    I did pip3 install ~/Downloads/python-igraph-0.7.1-1.tar.gz



## Test run
➜  hiconet python3 hiconet/HiCoNet.py hiconet/datasets/SDY80

    BTM conversion is automatically on for transcriptomics data (checking if gene number > 3000)

    To convert transcriptomics to BTM activity scores:
    ➜  Code cd hiconet/hiconet/btm
    ➜  btm git:(master) ✗ ls
    __init__.py         btm_example_data.py btm_tool.py
    ➜  btm git:(master) ✗ python
    >>> from btm_tool import genetable_to_activityscores
    >>> 
    >>> genetable_to_activityscores('/Users/sli/Desktop/SDY522/immport_results_SDY522_SDY522_Other_LAIV_Expression_Matrices.tsv', 
        '/Users/sli/Desktop/SDY522/BTMactivity_SDY522.txt')
    >>> 


## Planning
Need set up a web server asap (JR?)

Need js based visualization of data. 
Current version has HiCoNet.export_json working.

Build Docker container.

Apply to ImmPort datasets.

Add time_difference feature (SL)



