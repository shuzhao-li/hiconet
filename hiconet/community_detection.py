"""Hieracrchical Community Network
for integration of multiple data types

This module contains the community detection methods, all of which use
inputDataTable.dataframe as input, i.e. dataframe based on a matrix of samples by features.
Return format???
community structure, log

Local communities will be determined by community_detection methods, 
and be attached to the society class

Methods of community detection to consider:
a) HCL and its LC-MS variant
b) Girvan-Newman method, from networkx
c) Leidenalg, https://github.com/vtraag/leidenalg



input dataframe has to be correctly read: values start at position [0,0].
    Extra rows/columns go to SampleDictionary/FeatureDictionary


"""

__all__ = ["hierachical_clustering", "hierachical_clustering_lcms",
           #"louvain_find_communities",
           "leiden_find_communities",
            # to add other methods
          ]

import numpy as np
from scipy.cluster.hierarchy import *
from scipy.spatial.distance import pdist

import scanpy as sc
import anndata


def hierachical_clustering(df, distanceCut = 2):
    """Finding communities by hierachical clustering then cutting dendrogram tree.
    
    Parameters
    ----------
    df: input dataframe (data matrix)
    distanceCut is used to cut dendrogram tree. This can be optimized.

    Returns
    -------
    Clus: list of membership
    ClusDict: dictionary of communities to members
    """

    # distance matrix
    # print (df.values[:2, 1:5])
    # Y = pdist(df.values[:, 1:], 'correlation')
    Y = pdist(df.values, 'correlation')
    print(df.shape, Y.shape)

    # linkage matrix
    Z = linkage(Y, method='ward')
    Clus = fcluster(Z, distanceCut, criterion='distance')

    print(Clus)    # This is cluster number for each row in df

    number_features, number_clusters = len(Clus), len(set(list(Clus)))
    print("number of features: ", number_features)
    print("number of communities: ", number_clusters)

    # Compile clusters
    ClusDict = {}
    for ii in range(number_features):
        # if ClusDict.has_key(Clus[ii]):
        if Clus[ii] in ClusDict:
            ClusDict[ Clus[ii] ].append(ii)
        else:
            ClusDict[ Clus[ii] ] = [ii]

    #print(ClusDict.items()[:3])    # This organizes cluster, members
    return Clus, ClusDict



def girvan_newton_spectral_clustering():
    """
    An implementation is in mummichog code, to be phased out.
    """
    pass


def hierachical_clustering_lcms(DataMatrix, FeatureAnnotation, distanceCut = 3):
    """Clustering of LC-MS data by considering retention time.
    This uses the trio data structure
    
    Metabolomics FeatureAnnotation must have columns 'mz', 'rtime'.
    
    porting into HiCoNet
    """

    # Clustering of metabolite features
    # distance matrix, this is [1 - (Pearson R)]
    metabo = DataMatrix
    YM = pdist(metabo.values[:, 1:], 'correlation')
    print(metabo.shape, YM.shape)

    # New method, weighting delta retention time into new distance matrix
    retention_time = [FeatureAnnotation.loc[x, 'rtime'] for x in metabo.index]         
    # Note this relies on correct reading table; index is mz_rt
    min_retention_time, max_retention_time = min(retention_time), max(retention_time)
    range_retention_time = max_retention_time - min_retention_time
    print("min_retention_time, max_retention_time", min_retention_time, max_retention_time)


    PearsonR = 1 - YM
    #print(PearsonR[:20])

    delta_RT = []
    for ii in range(metabo.shape[0]):
        for jj in range(ii+1, metabo.shape[0]):
            delta_RT.append(abs(retention_time[ii] - retention_time[jj]))

    print("Vector delta_RT len: ", len(delta_RT))

    #
    # weighting function
    # distance = 1 - (1 - delta_RT/range_retention_time)*PearsonR
    #

    YM_new = 1 - (1- np.array(delta_RT)/range_retention_time)*PearsonR
    ZM = linkage(YM_new, method='ward')
    metClus = fcluster(ZM, distanceCut, criterion='distance')
    #print(metClus[:10])

    number_features, number_clusters = len(metClus), len(set(list(metClus)))
    print("number of features: ", number_features)
    print("number of communities: ", number_clusters)

    # Compile clusters
    metClusDict = {}
    for ii in range(number_features):
        # if metClusDict.has_key(metClus[ii]):
        if metClus[ii] in metClusDict:
            metClusDict[ metClus[ii] ].append(ii)
        else:
            metClusDict[ metClus[ii] ] = [ii]

    return metClus, metClusDict


def leiden_find_communities(df, method='modularity'):
    """
    
    https://github.com/vtraag/leidenalg
    This finds the optimal partition using the Leiden algorithm [1], 
    which is an extension of the Louvain algorithm [2] for a number of different methods. 
    leidenalg.find_partition(G, leidenalg.ModularityVertexPartition);
        
    Parameters
    ----------
    df: input dataframe (data matrix)
    method: default measure modularity

    Returns
    -------
    Clus: list of membership
    ClusDict: dictionary of communities to members
    """
    #
    scdm = anndata.AnnData(df)
    sc.pp.neighbors(scdm, use_rep='X')
    sc.tl.leiden(scdm)
    Clus = list(scdm.obs['leiden'].values)
    #
    # using old code below. To be optimized
    #
    number_features, number_clusters = len(Clus), len(set(list(Clus)))
    print("number of features: ", number_features)
    print("number of communities: ", number_clusters)

    # Compile clusters
    ClusDict = {}
    for ii in range(number_features):
        # if ClusDict.has_key(Clus[ii]):
        if Clus[ii] in ClusDict:
            ClusDict[ Clus[ii] ].append(ii)
        else:
            ClusDict[ Clus[ii] ] = [ii]

    #print(ClusDict.items()[:3])    # This organizes cluster, members
    return Clus, ClusDict


    
def louvain_find_communities(df, method='modularity'):
    pass