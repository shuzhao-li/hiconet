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


    Graph: as in loom (loompy.org).
    Community: a group of features within a society that share a similar pattern.
    e.g. one way to identify communities is to do Louvain graph-clustering method on the graph.


    Reusing
    -------
    The data structure is aligned with anndata (https://github.com/theislab/anndata) but transposed.
    For loom format, one loom file = a Society; as meta data can differ for different data types.
    Same goes for anndata
