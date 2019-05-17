'''
Created on 12.07.2017

@author: marisakoe
'''

import codecs
from collections import defaultdict

###################data method####################################
def get_data_dyen():
    '''
    reads the data from willems et al and creates the data matrix and the distance matrix for the hybridization network computation
    '''
    
#     f = open("input/dyen_database.csv")
#     rawlines = f.readlines()
#     f.close()
    
    f = codecs.open("input/dyen_database.csv","r","utf-8")
    rawlines=f.readlines()
    f.close()
    
    ##dictionary with sets of languages for meaning m
    L_m = defaultdict(list)
    ##dictionar with sets of cognate sets for meaning m
    C_m = defaultdict(list)
    ##dict with key=meaning, value=dict with key=lang value=[cc]
    langs_cc = defaultdict(lambda: defaultdict(list))
    ##dict with key=meaning, value = dict with key=lange value=[words]
    langs_words = defaultdict(lambda: defaultdict(list))
    ##go through the lines and read the data
    for l in rawlines[1:]:
        line = l.split(",")
        concept = line[0]
        cc = line[1]
        
        lang = line[2]
        #langNr = line[3]
        word = line[4].split("\n")[0].split("\r")[0] 
        
        concat_word = [word,cc]
        
        if lang not in L_m[concept]:
            L_m[concept].append(lang)
        if cc not in C_m[concept]:
            C_m[concept].append(cc)
        langs_cc[concept][lang].append(cc)
        langs_words[concept][lang].append(concat_word)
        
    
    
    return L_m, C_m, langs_cc, langs_words



##########writing methods##################################

def write_dm(distance_matrix, folderFileName):
    '''
    writes the distance matrix in a file
    :param distance_matrix:the distance matrix
    '''
    
    ##get all the languages in the sample
    keys = distance_matrix.keys()
    
    fout = codecs.open("output/dyen/"+folderFileName+".phy","wb","utf-8")
    ##writes the number of the languages at the top of the file
    fout.write(str(len(keys))+"\n")
    ##writes the distances in the matrix in the phylib format
    for lang, valdict in distance_matrix.items():
        row = lang.ljust(20)+"\t"
        
        for lang2 in valdict:
            row=row+str(distance_matrix[lang][lang2])+" "
            
        fout.write(row+"\n")
        
    fout.close()

if __name__ == '__main__':
    pass