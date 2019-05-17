'''
Created on 12.07.2017

@author: marisa
'''

import codecs, itertools, distance
from collections import defaultdict

from helper_methods import get_data_nelex, read_nelexPMI, write_dm

def nelex_compute_dm_levenshtein():
    '''
    computes the distance matrices for each concepts and saves them in a folder.
    the formular is according to Willems et al 
    d_m(l1,l2) = sum(dc(l1,l2))/n_l1l2 
    the distance matrices are used to compute concept trees wiht R
    '''
    L_m, C_m, m_lang_cc, m_lang_words = read_nelexPMI()
    
    ##for each meanidn in the dict, get the list of languages
    for meaning, list_langs in L_m.items():
        #print meaning
        #if meaning == "Berg::N":
        ##create the distance matrix for this meaning and set the language pairs 0.0
        distance_matrix = defaultdict(lambda: defaultdict(float))
        combineList = list(itertools.product(list_langs,list_langs))
        for pair in combineList:
            distance_matrix[pair[0]][pair[1]] = 0.0
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
         
    
        ##get the combinations for the languages
        list_lang_pairs = itertools.combinations(list_langs, 2)
        ##for each language pair
        for pair in list_lang_pairs:
            ##get the list of cognate classes for each language
            l1_cc = data_matrix[pair[0]]
            l2_cc = data_matrix[pair[1]]
            ##get the list of words plus cognate classes for each language
            l1_words = m_lang_words[meaning][pair[0]]
            l2_words = m_lang_words[meaning][pair[1]]
            ##calculate the distance between the languages
            lang_dist = calculate_dc(l1_cc, l2_cc, l1_words, l2_words)
            ##fill the distance matrix with the values
            distance_matrix[pair[0]][pair[1]] = lang_dist
            distance_matrix[pair[1]][pair[0]] = lang_dist
        ##choose a file name
        fileName = "dstmtxLev/distanceMatrixLev_"+meaning
        ##write the distance matrix to the output folder
        write_dm(distance_matrix, fileName)
        
        
###################computation methods#############################

def calculate_dc(l1_cc,l2_cc, l1_words, l2_words):
    '''
    Calculate the distance between two languages in one concept.
    sum(dc(l1,l2))/n_l1,l2
    where sum(dc(l1,l2)) is the sum over all cogante in the set for this concept and dc is the distance for the cognat between l1 and l2
    dc = 0 if neither word forms were present in c
    dc = 1 if either word l1 or word l2 are present in c
    dc = min(levensthein) if both have a word for this cognate class
    :param l1_cc:dict with key=cognate class value=binary 1 for presence 0 for absence
    :param l2_cc:dict with key=cognate class value=binary 1 for presence 0 for absence
    :param l1_words:words for language 1
    :param l2_words:words for language 2
    '''
    #print l1_cc, l2_cc, l1_words, l2_words
    
    
    
    list_dist = []
    n_l1l2 = 0.0
    for cc, l1_val in l1_cc.items():
        ##if both values are 0 (absence of cogante class for the language) the distance is 0
        if l1_val==0 and l2_cc[cc]==0:
            
            list_dist.append(0.0)
        ##if both values are 1 (presence of cognate class for the languages) the distance is computed with the normalized Levenshtein distance
        elif l1_val==1 and l2_cc[cc]==1:
            ##get the list of words for each language
            l1_cc_words = []
            l2_cc_words=[]
            ##filter out the words in the same cogante class for the computation
            for word in l1_words:
                if word[1] == cc:
                    l1_cc_words.append(word[0])
            for word2 in l2_words:
                if word2[1]==cc:
                    l2_cc_words.append(word2[0])
            
            ##compute the levensthtein distance
            dist =compute_levenshtein(l1_cc_words, l2_cc_words)
            ##append it to the list (because we are taking the sum)
            list_dist.append(float(dist))
            ##compute n_l1l2
            n_l1l2 += 1
        ##if either l1 or l2 are in the cognate data the distance is 1.0
        else:
            list_dist.append(1.0)
            n_l1l2 += 1
            
    ##compute the formular and return the distance
    end_dist = sum(list_dist)/n_l1l2
    return end_dist

def compute_levenshtein(l1_words, l2_words):
    '''
    Compute the normalized Levenshtein distance for all words in the same cognate class.
    It can take synonyms into account, where we will take the minimal value between the lists.
    :param l1_words:list of words for l1
    :param l2_words:list of words for l2
    :return min(lev_dist): the minimum value of the list of distances (only important for synonyms)
    '''
    lev_dist = []
    combineList = list(itertools.product(l1_words,l2_words))
    for pair in combineList:
        ##normalized levenshtein distance the longer alignment is the factor (method=1 takes the shorter alignment as factor)
        dist = distance.nlevenshtein(pair[0], pair[1], method=2)
        lev_dist.append(dist)
    ##take the minimal value of the list if there are synonyms, otherwise there is only one value
    return min(lev_dist)




if __name__ == '__main__':
    pass