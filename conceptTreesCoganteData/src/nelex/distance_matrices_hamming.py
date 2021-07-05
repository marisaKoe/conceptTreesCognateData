'''
Created on 12.07.2017

@author: marisakoe

NOTE: this is independent of the other distance measures. This was just for testing


'''
import codecs, itertools, scipy.spatial
from collections import defaultdict

from helper_methods import get_data_nelex, write_dm

def nelex_compute_dm_hamming():
    '''
    reads the data from willems et al and creates the data matrix and the distance matrix 
    '''
    
    L_m, C_m, m_lang_cc, m_lang_words = get_data_nelex()
    
    for meaning, list_langs in L_m.items():
        ##create the distance matrix for this meaning and set the language pairs 0.0

        ##create a small datamatrix for the computation and set everything default 0
        data_matrix = defaultdict(lambda: defaultdict(int))
        combineList = list(itertools.product(list_langs,C_m[meaning]))
        for pair in combineList:
            data_matrix[pair[0]][pair[1]] = 0
        
        ##for each lang, dic in the data_matrix fill the matrix
        for lang, cc_dict in data_matrix.items():
            ##get the list of cognate classes for the concept and the language
            lang_cc = m_lang_cc[meaning][lang]
            ##get all the cogante classes for this meaning
            cc_keys = cc_dict.keys()
            ##for each cogante class for the language, check if it is in keys and set to 1
            for cc in lang_cc:
                if cc in cc_keys:
                    data_matrix[lang][cc] = 1

    

    
        ##compute the distance matrix for willems et al with two versions
        fileName1 = "dstmtxHam/hamming/willemsMatrix_"+meaning
        fileName2 = "dstmtxHam/normHamming/willemsNormMatrix_"+meaning
        version = "hamming"
        version_norm = "norm_hamming"
        compute_dm(data_matrix, list_langs, fileName1, version)
        compute_dm(data_matrix, list_langs, fileName2, version_norm)
    
#############matrix methods######################################


                    
        
def compute_dm(data_matrix, langs, fileName, version):
    '''
    compute a distance matrix out of a data matrix with binary characters using the Hamming distance.
    :param data_matrix: the data matrix dict with key=langs val=dict with key=cc value= 0 for absence 1 for presence
    :param langs: list with all unique languages
    :param fileName: the name for the output file
    '''
    ##create the default distance matrix
    distance_matrix = defaultdict(lambda: defaultdict(float))
    combineList = list(itertools.product(langs,langs))
    for pair in combineList:
        distance_matrix[pair[0]][pair[1]] = 0.0
        
    ##get the languages in the data_matrix and the values as arrays to comput the hamming distance between two languages
    for pair in combineList:
        if not pair[0] == pair[1]:
            val1 = data_matrix[pair[0]].values()
            val2 = data_matrix[pair[1]].values()
            if version == "hamming":
                ##hamming distance
                dis = hamdist(val1, val2)
                distance_matrix[pair[0]][pair[1]] = dis
                
            elif version == "norm_hamming": 
                ##normalized hamming distance
                dis = scipy.spatial.distance.hamming(val1, val2)
                distance_matrix[pair[0]][pair[1]] = dis
    
    ##write the matrix
    write_dm(distance_matrix, fileName)
    #return distance_matrix
    

######################helper methods########################


def hamdist(str1, str2):        
    '''
    Compute the Hamming distance between two strings
    :param str1:
    :param str2:
    '''
    diffs = 0.0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs+=1
    
    return diffs


            

if __name__ == '__main__':
    pass
    
    
    