'''
Created on 26.07.2017

@author: marisakoe
'''

import codecs, itertools, distance
from collections import defaultdict
from numpy import *

from helper_methods import read_nelexPMI, write_dm

def nelex_compute_dm_nw():
    '''
    computes the distance matrices for each concepts and saves them in a folder.
    the formular is according to Willems et al 
    d_m(l1,l2) = sum(dc(l1,l2))/n_l1l2 
    the distance matrices are used to compute concept trees wiht R
    '''
    L_m, C_m, m_lang_cc, m_lang_words = read_nelexPMI()
    
    
    ##for each meaning in the dict, get the list of languages
    for meaning, list_langs in L_m.items():
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
        fileName = "dstmtxNW/distanceMatrixNW_"+meaning
        ##write the distance matrix to the output folder
        write_dm(distance_matrix, fileName)



####################computation methods#############################

gp1 = -2.49302792222
gp2 = -1.70573165621

lodict={}


def calculate_dc(l1_cc,l2_cc, l1_words, l2_words):
    '''
    Calculate the distance between two languages in one concept.
    sum(dc(l1,l2))/n_l1,l2
    where sum(dc(l1,l2)) is the sum over all cognates in the set for this concept and dc is the distance for the cognate between l1 and l2
    dc = 0 if neither word forms were present in c
    dc = 1 if either word l1 or word l2 are present in c
    dc = min(levensthein) if both have a word for this cognate class
    :param l1_cc:dict with key=cognate class value=binary 1 for presence 0 for absence
    :param l2_cc:dict with key=cognate class value=binary 1 for presence 0 for absence
    :param l1_words:words for language 1
    :param l2_words:words for language 2
    '''

    #reads the file with the ASJP sounds present in nelex
    f = open('input/sounds.txt')
    sounds = array([x.strip() for x in f.readlines()])
    f.close()
    
    #reads the log odds scores for the whole word languages
    f = open('input/pmi_matrix_nelex.txt','r')
    l = f.readlines()
    f.close()
    logOdds = array([x.strip().split() for x in l],double)
    
    #assigns the log odds to the sound pairs to create the lodict
    for i in xrange(len(sounds)):#Initiate sound dictionary
        for j in xrange(len(sounds)):
            lodict[sounds[i],sounds[j]] = logOdds[i,j]
    
    
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

            #compute the needleman wunsch score and alignment for the words of the languages
            dist = nw(l1_cc_words, l2_cc_words, lodict, gp1, gp2)
            
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


    
    
    
    
def nw(l1_words,l2_words,lodict,gp1,gp2):
    """
    Needleman-Wunsch algorithm for pairwise string alignment
    with affine gap penalties.
    'lodict' must be a dictionary with all symbol pairs as keys
    and match scores as values.
    gp1 and gp2 are gap penalties for opening/extending a gap.
    Returns the alignment score and one optimal alignment.
    """
    nw_dist=[]
    combineList = list(itertools.product(l1_words,l2_words))
    for pair in combineList:
    
        x = pair[0]
        y= pair[1]
        #length of the words
        n,m = len(x),len(y)
        dp = zeros((n+1,m+1))
        pointers = zeros((n+1,m+1),int)
        for i in xrange(1,n+1):
            dp[i,0] = dp[i-1,0]+(gp2 if i>1 else gp1)
            pointers[i,0]=1
        for j in xrange(1,m+1):
            dp[0,j] = dp[0,j-1]+(gp2 if j>1 else gp1)
            pointers[0,j]=2
        for i in xrange(1,n+1):
            for j in xrange(1,m+1):
                match = dp[i-1,j-1]+lodict[x[i-1],y[j-1]]
                insert = dp[i-1,j]+(gp2 if pointers[i-1,j]==1 else gp1)
                delet = dp[i,j-1]+(gp2 if pointers[i,j-1]==2 else gp1)
                dp[i,j] = max([match,insert,delet])
                pointers[i,j] = argmax([match,insert,delet])
        alg = []
        i,j = n,m
        while(i>0 or j>0):
            pt = pointers[i,j]
            if pt==0:
                i-=1
                j-=1
                alg = [[x[i],y[j]]]+alg
            if pt==1:
                i-=1
                alg = [[x[i],'-']]+alg
            if pt==2:
                j-=1
                alg = [['-',y[j]]]+alg
        nw_dist.append(dp[-1,-1])
    
    ##take the minimal value of the list if there are synonyms, otherwise there is only one value
    #return max(nw_dist)
    
    #wrd_score = 1.0-((2.0*nw(word_pair[0],word_pair[1]))/(nw(word_pair[0],word_pair[0])+nw(word_pair[1],word_pair[1])))
    
    
    ##sigmoid function zur Berechnung der Distanz
    ##Vorteil Sigmoid ist eine S-Kurve. Alle Werte unter 0.5 signalisieren eine negative Aehnlichkeit (-> nicht aehnlich), all Werte ueber 0.5 
    ## signalisieren eine positive Aehnlichkeit (-> aehnlich). Dadurch kann der threshold fuer andere funktionen (lingpy) immer 0.5 sein.
    dist = 1.0 - (1.0/(1 + exp(-max(nw_dist))))
    return dist
    #return dp[-1,-1],alg,array([''.join(x) for x in array(alg).T])





if __name__ == '__main__':
    pass