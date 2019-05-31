# Hieracrchical Community Network (HiCoNet)
for integration of multiple data types collected from a common group of subjects.

An -omics dataset contains a lot of redundancy, and features of similar quantitative patterns can be considered as communities. Common methods of feature-level integration may exacerbate the problem of redundancy, as the combination space gets large and complex. HiCoNet detects communities within each data type, then tests the association between communities across data types. This "hierarchical community network" thus provides a reasonable model of the organizational structure of measured biology.

We started this approach in a VZV vaccine study (Li et al, 2017, https://doi.org/10.1016/j.cell.2017.04.026), and the method was further developed through other scientific projects.

## The 3-file-society Data Strucutre
Each data type often has its own idiosyncrasy. To be able to automate sophisticated analysis, a most generic format of common denominator is desired.
We use three files to describe one data type, DataMatrix, FeatureAnnotation and ObservationAnnotation.
The DataMatrix file uses a single row for observation IDs and a single column for feature IDs. This mandates unique identifier per feature per observation, and separate meta data from the DataMatrix. Because the annotations on feature or observation can be heterogenuous, but should not affect the format of DataMatrix. 

## More on Terminology:
Study: an administrative unit that include one or more projects. Same as ImmPort "Study" (https://www.immport.org/resources/documentation). 

Project: a collection of data of one or more types (a dataset). For multiple data types, common samples/subjects are expected, as HiCoNet deals with the `N-integration` problem.
    This is the unit HiCoNet works on - HiCoNet integrates DataMatrices within a DataSet
    A DataSet should have at least one Society of data.

Society: one data type, defined by a set of DataMatrix, FeatureAnnotation (optional) and ObservationAnnotation (optional). 
    This 3-file design is similar to anndata (https://github.com/theislab/anndata) but data are transposed. Meta data can differ for different data types.

DataMatrix: a data matrix of [continuous] values that represent a biological state or concentration, of the same data type.
    E.g. transcriptomics (array intensity or transcript counts), metabolomics (peak intensity/area) or microbial OTU counts.
    This can include different time points or treatments. This is the unit community detection is based on.

ObservationAnnotation: an observation is an experimental measurement of a biological sample. 
    A sample may have measurement replicates. Description of biological samples should be in ObservationAnnotation, which can support inferring the study design (e.g. treatment, time points). For ImmPort data, the MySQL table `biosample` can serve as ObservationAnnotation. Time points and treatment are key annotation variables in many studies. 

FeatureAnnotation: meta data on features. 
    This can be as simple as gene annotation, which can even be optional. But a feature may carry a defition of multiple parameters. E.g. a metabolite feaure may have m/z, retention time and collision cross section, and these parameters may be used for certain algorithms.

Graph: a graph/network for relationships in the data (e.g. used in loom format, loompy.org). The current version of HiCoNet does not store this, but will consider it for future versions.

Community: a group of features within a society that share a similar pattern.

## Requires
    'PyYAML'
    'numpy',
    'scipy',
    'pandas',
    'sklearn',
    'leidenalg',
    'scanpy',
    'igraph',
    'fuzzywuzzy',

    Note: python-igraph requires the C library igraph. The installation on Mac OS may be tricky:  
    https://stackoverflow.com/questions/45667147/install-python-igraph-on-mac
    I did pip3 install ~/Downloads/python-igraph-0.7.1-1.tar.gz
    For a Docker or new install, both igraph and python-igraph are needed.

## Use
    This software package is available via PyPI (Python Package Index) and GitHub.
    Test datasets are included. E.g. to run test:
    python3 -m hiconet.HiCoNet hiconet/datasets/SDY80

    There are related but separate projects of hiconet-server and hiconet-explorer.
    