"""
Hieracrchical Community Network 
for integration of multiple data types

This module contains the function to build 
PLS2 association networks 
"""

from random import sample as permutation
import numpy as np
from sklearn.cross_decomposition import PLSRegression

from .input_functions import NUM_PERMUTATION


class pairNetwork:
    """
    This is a unit to perform PLS regression on two data types (societies).
    The societies can be sliced by attributes; 
    communities were pre-computed at the initiation of Society class.
    Permutation is performed within this class.
    Result is stored as a network of associated communities.


    observation_list_society1, observation_list_society2, expected in association_dict,
    which should have done sample/subj matching 

    Permutation is done here to compute p-values for association R.
    Resampling is done on the whole society.

    """
    def __init__(self, association_dict, society1, society2):
        # 
        self.dict = association_dict
        self.name = association_dict['name'] + '_' + str(association_dict['timepoint1']) + '_' + str(association_dict['timepoint2'])
        self.SampleNumber = len( association_dict['observation_list_society1'])
        
        # this runs PLS2 and get p-values via permutation test
        self.PLS_2datatypes(association_dict, society1, society2)
        # network_edges = [( g, m, PLSscore, p-value ), ...], g and m as community IDs from society1/2
        self.network_edges = self.pls_network_edges
     

    def PLS_2datatypes(self, association_dict, society1, society2, ):
        '''
        Subjs, samples matched per instruction from association_dict
    
        '''
        label1, label2 = society1.name, society2.name
        # g and m here are legacy names, without specific meaning here
        gCommunities, mCommunities = society1.Communities, society2.Communities
        gSizes = [len(x) for x in gCommunities.values() if len(x) > 2]
        mSizes = [len(x) for x in mCommunities.values() if len(x) > 2]
        
        gArray, mArray = list(society1.DataMatrix.values.flatten()), \
                         list(society2.DataMatrix.values.flatten())
        gDF, mDF = society1.DataMatrix[association_dict['observation_list_society1']], \
                   society2.DataMatrix[association_dict['observation_list_society2']]
        
        # get PLS2 scores
        pls_scores = self.get_pls_scores_real(gCommunities, mCommunities, gDF, mDF)
        # Do permutation
        permutation_scores = self.get_pls_scores_permutation(gArray, mArray, gSizes, mSizes)
        
        #print("pls_scores", pls_scores)
        #print("permutation_scores", permutation_scores)
        
        # Collect network edges [( node1, node2, PLSscore, p-value ), ...]
        self.pls_network_edges = self.get_p_values(pls_scores, permutation_scores)
        
        
    def get_p_values(self, pls_scores, permutation_scores):
        '''
        Compute p-value by polyfit the permutation_scores
        
        Parameters
        ----------
        pls_scores, permutation_scores
        
        Returns
        -------
        p-value list, [( g, m, PLSscore, p-value ), ...], g and m as community IDs from society1/2
        '''
        # polyfit permutation scores
        # look up p-values for PLS scores via a fitted polynomial function. Rank based method is too slow. 
        # More interested in the top 50%, and should down-sample in future version
        plslistlength = len(permutation_scores)
        # not limited to positive score here, to be more robust, i.e. tolerant to crapy data
        permutation_scores_positive = permutation_scores    #[x for x in permutation_scores if x > 0]
        N_positive = len(permutation_scores_positive)
        N_to_use = int(N_positive/2)
        permutation_scores_positive.sort(reverse=True)
    
        # fit a function log10(p-value) = f(score)
        fp = np.polyfit( permutation_scores_positive[:N_to_use], 
                         np.log10(np.arange(1,N_to_use+1)/float(plslistlength)), 
                         3)
    
        def fp_predict(x, fp):
            return fp[0]*x**3 + fp[1]*x**2 + fp[2]*x + fp[3]
    
        newlist = []
        for x in pls_scores:
            # ( g, m, PLSscore ), PLSscore is x[2]
            # force p <= 1
            newlist.append( (x[0], x[1], x[2], 
                             min(1, 10**fp_predict(x[2], fp))) )
        # sort by PLSscore decending
        def sort2(val): return val[2]
        newlist.sort(key = sort2, reverse = True)
        return newlist
    
    
    def write_p_values(self, s):
        #
        # Not functional, not used for now
        #
        # datatype_1\tcommunity_number\tdatatype_2\tcommunity_number\tPLS_score\tp-value\n
        #society1, observation_list_society1, society2, observation_list_society2, network_name

        outFile = label1 + '_' + label2 + '__PLSscore.txt'
        print(outFile)
        s = '\t'.join( [label1, label2, 'PLSscore', 'log10(p-value)'] ) + '\n'
        for x in pls_scores:
            # ( g, m, PLSscore ), PLSscore is x[2]
            if x[2] > 0:
                log10p = fp_predict(x[2], fp)
                line = '\t'.join([str(ii) for ii in x] + [str(log10p)])
                if log10p < -2:  print( line )
                s += line + '\n'
    
        with open(outFile, 'w') as file:
            file.write(s)
            
        
        
    def plot_polyfit_scores(self, plslistlength, permutation_scores_positive):
        #
        # Not used for now
        #
        # plot permutation polyfit
        plt.figure()
        step = int(N_to_use/100)
        sparseY = np.arange(1, N_to_use, step)
        sparseX = [permutation_scores_positive[ii] for ii in sparseY]
        plt.plot(sparseX, np.log10(sparseY/float(plslistlength)), 'bo')
        plt.plot(sparseX, [fp_predict(x, fp) for x in sparseX], 'r-')
        plt.savefig("polyfit_" + outFile.replace(".txt", ".png"))
        
        
        

    def get_pls_scores_real(self, gCommunities, mCommunities, gDF, mDF):
        '''
        Compute PLS2 scores for all pairwise communities from two societies.
        
        Parameters
        ----------
        gCommunities, mCommunities, gDF, mDF
        Communities and DataMatrix from society_1, society_2
        
        Returns
        -------
        pls_scores list as [( g, m, PLSscore ), ...]
        '''
        PLS = PLSRegression(n_components=3)
        pls_scores = []
    
        for g in gCommunities.keys():
            if len(gCommunities[g]) >= 3:
                #print(g,)
                for m in mCommunities.keys():
                    if len(mCommunities[m]) >= 3:
                        # Getting corresponding rows from btm and metabo. 
                        matrix1, matrix2 = gDF.values[ gCommunities[g], : ], mDF.values[ mCommunities[m], : ]
                        matrix1, matrix2 = np.transpose(matrix1), np.transpose(matrix2)
                        
                        print("input matrices ", matrix1.shape, matrix2.shape)
                        # PLS regression
                        if matrix1.shape[1] > matrix2.shape[1]:
                            PLS.fit(matrix1, matrix2)
                            PLSscore = PLS.score(matrix1, matrix2) 
                        else:
                            PLS.fit(matrix2, matrix1)
                            PLSscore = PLS.score(matrix2, matrix1)
    
                        pls_scores.append(( g, m, PLSscore ))
    
        return pls_scores


    def get_pls_scores_permutation(self, gArray, mArray, gSizes, mSizes, 
                                   numPermutation=NUM_PERMUTATION):
        '''
        g and m from legacy code, no particular meaning
        
        ???
        Permutation will be done within each data slice,
        due to different data characteristics in time points or delta, or etc.
        
        '''
        SampleNumber = self.SampleNumber
        PLS = PLSRegression(n_components=3)
        scores = []
        for jj in range(numPermutation):
            if str(jj)[-1] == '0': print ("            Permutation --- %d" %jj)
    
            for g in gSizes:
                matrix1 = []
                for ii in range(g):
                    matrix1.append(permutation(gArray, SampleNumber))
                matrix1 = np.array(matrix1).T
                for m in mSizes:
                    matrix2 = []
                    for ii in range(m):
                        matrix2.append(permutation(mArray, SampleNumber))
                    matrix2 = np.array(matrix2).T
                    #
                    if matrix1.shape[1] > matrix2.shape[1]:
                        PLS.fit(matrix1, matrix2)
                        PLSscore = PLS.score(matrix1, matrix2)
                    else:
                        PLS.fit(matrix2, matrix1)
                        PLSscore = PLS.score(matrix2, matrix1)
    
                    scores.append( PLSscore )
    
        return scores
    
