'''
Created on 12.07.2017

@author: marisa
'''
import codecs
from collections import defaultdict

def get_data_nelex():
    '''
    reads the file with the cognates
    :return L_m: dictionary with sets of languages for meaning m key=concept value=[langs]
    :reutrn C_m: dictionary with sets of cognate sets for meaning m key=concept value=[cognagte sets]
    :return langs_cc: dict with key=concept value=dict with key=lang value=[cc]
    :return langs_words: dict with key=concept value=dict with key=lang value=[words]
    '''
    
    ##read the nelex data
    f = codecs.open("input/inferred-cognates-northeuralex0.9-single-thresh0.55-etym.tsv","r","utf-8")
    rawlines = f.readlines()
    f.close()
    
    ##dictionary with sets of languages for meaning m
    L_m = defaultdict(list)
    ##dictionary with sets of cognate sets for meaning m
    C_m = defaultdict(list)
    ##dict with key=meaning, value=dict with key=lang value=[cc]
    langs_cc = defaultdict(lambda: defaultdict(list))
    ##dict with key=meaning, value = dict with key=lange value=[words]
    langs_words = defaultdict(lambda: defaultdict(list))
    
    
    ## 0=concept, 1=iso, 2=originalform, 3=ipa, 4=cc, 5=loan
    for line in rawlines:
        line = line.split("\t")
        
        ##concept
        concept = line[0]

        if " " in concept:
            #print concept
            concept = concept.replace(" ","")
            #print concept
        ##iso lang
        lang = line[1]
        ##originalform
        #orig_word = line[2]
        ##ipa form
        ipa_wrd = line[3]
        ##cognate class
        cc = line[4]
        ##loan info
        #loan = line[-1]
        
        concat_word = [ipa_wrd,cc]
        
        if lang not in L_m[concept]:
            L_m[concept].append(lang)
        if cc not in C_m[concept]:
            C_m[concept].append(cc)
        langs_cc[concept][lang].append(cc)
        langs_words[concept][lang].append(concat_word)
        
    
    
    return L_m, C_m, langs_cc, langs_words

def read_nelexPMI():
    '''
    reads the nelex datafile, which can be downloaded from the website with the ASJP transcriptions
    read the file created with Tarakas program, since we need cognates
    :return L_m: dictionary with sets of languages for meaning m key=concept value=[langs]
    :reutrn C_m: dictionary with sets of cognate sets for meaning m key=concept value=[cognagte sets]
    :return langs_cc: dict with key=concept value=dict with key=lang value=[cc]
    :return langs_words: dict with key=concept value=dict with key=lang value=[words]
    '''
    
    
    f1 = open("input/nelexAsjp.cognates",'r')
    raw_data = f1.readlines()
    f1.close()
    
    ##dictionary with sets of languages for meaning m
    L_m = defaultdict(list)
    ##dictionary with sets of cognate sets for meaning m
    C_m = defaultdict(list)
    ##dict with key=meaning, value=dict with key=lang value=[cc]
    langs_cc = defaultdict(lambda: defaultdict(list))
    ##dict with key=meaning, value = dict with key=lange value=[words]
    langs_words = defaultdict(lambda: defaultdict(list))


    for line in raw_data[1:]:
        line = line.strip()
        #split it by tab
        arr = line.split("\t")
        concept = arr[0]
        lang=arr[1]
        asjp=arr[2]
        cc = arr[-1]
        
        
        #word transcription
        asjp_word = arr[2].split(",")[0]
        asjp_word = asjp_word.replace(" ", "")
        #tokenized_word = ipa2tokens(asjp_word)
        #asjp_word = "".join(tokens2class(tokenized_word, 'asjp'))
        asjp_word = asjp_word.replace("%","")
        asjp_word = asjp_word.replace("~","")
        asjp_word = asjp_word.replace("*","")
        asjp_word = asjp_word.replace("$","")
        asjp_word = asjp_word.replace("\"","")
        asjp_word = asjp_word.replace(""" " ""","")
        asjp_word = asjp_word.replace("K","k")
        if len(asjp_word) < 1:
            continue
        if "K" in asjp_word:
            print line
            
        concat_word = [asjp,cc]
        
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
    
    fout = codecs.open("output/nelex/"+folderFileName+".phy","wb","utf-8")
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
    
    get_data_nelex()