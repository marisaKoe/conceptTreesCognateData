'''
Created on 12.07.2017

@author: marisakoe
'''

from dyen import dyen_compute_dm_levenshtein, dyen_compute_dm_hamming
from nelex import nelex_compute_dm_levenshtein, nelex_compute_dm_hamming, nelex_compute_dm_nw

import glob
import rpy2.robjects as r
from rpy2.robjects.packages import importr

#import the basis of R
base = importr("base")
utils = importr("utils")
stats = importr("stats")
#imports ape, required for phangorn
ape = importr("ape")
#imports phangorn
phangorn = importr("phangorn")

def reconstruct_trees_phy(mtx_folder, outputFolder):
    '''
    Reconstructs Neighbor Joining trees for each concept
    The trees are reconstructed using R and the packages ape and phangorn
    '''

    list_matrices = glob.glob(mtx_folder)
    print list_matrices
    for f in list_matrices:
        concept = f.split("/")[-1].split(".")[0].split("_")[-1]
        print concept
        t = utils.read_table(f, skip=1, row_names=1)
        mx = base.as_matrix(t)
        dm = stats.as_dist(mx)
        #nj trees
        tree = ape.nj(dm)
        ape.write_tree(tree, file=outputFolder+"nj/"+concept+"+njTree.nwk")
        #fastme trees
        tree1 = ape.fastme_bal(dm, nni=True, spr=True, tbr=False)
        ape.write_tree(tree1, file=outputFolder+"fastme/"+concept+"+fastmeTree.nwk")

if __name__ == '__main__':
    
    ##compute the matrices for the dyen data
    #dyen_compute_dm_levenshtein()
    #dyen_compute_dm_hamming()
    ##compute the trees
    ##levenshtein
#     lev_mtx = "output/dyen/dstmtxLev/*.phy"
#     outputLev = "output/dyen/conceptTrees/Levenshtein/"
#     reconstruct_trees_phy(lev_mtx, outputLev)
    ##hamming
#     ham_mtx = "output/dyen/dstmtxHam/hamming/*.phy"
#     outputHam = "output/dyen/conceptTrees/Hamming/"
#     reconstruct_trees_phy(ham_mtx, outputHam)
#     #list_matrices = glob.glob("../output/dyen/dstmtxHam/hamming/*.phy")
#     ##normalized hamming
#     normham_mtx = "output/dyen/dstmtxHam/normHamming/*.phy"
#     outputNormHam = "output/dyen/conceptTrees/NormHamming/"
#     reconstruct_trees_phy(normham_mtx, outputNormHam)
    
    
    
    
    ##compute matrices for nelex
    
    ##nw with pmi, sounds and cognates from Tarakas method
    
    nelex_compute_dm_nw()
    nw_mtx = "output/nelex/dstmtxNW/*.phy"
    outputNW = "output/nelex/conceptTrees/NW/"
    reconstruct_trees_phy(nw_mtx, outputNW)

#     ##levenshtein
#     nelex_compute_dm_levenshtein()
#     lev_mtx = "output/nelex/dstmtxLev/*.phy"
#     outputLev = "output/nelex/conceptTrees/Levenshtein/"
#     reconstruct_trees_phy(lev_mtx, outputLev)
#     ##hamming
#     nelex_compute_dm_hamming()
#     ham_mtx = "output/nelex/dstmtxHam/hamming/*.phy"
#     outputHam = "output/nelex/conceptTrees/Hamming/"
#     reconstruct_trees_phy(ham_mtx, outputHam)
#     ##normalized hamming
#     normham_mtx = "output/nelex/dstmtxHam/normHamming/*.phy"
#     outputNormHam = "output/nelex/conceptTrees/NormHamming/"
#     reconstruct_trees_phy(normham_mtx, outputNormHam)
    
    
    
    
    
    