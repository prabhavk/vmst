# Implementation of the minimum leaves vmst algorithm 

from vmst_library import ConstructMinLeavesVMST
# Edge weights = pairwise distances [upper/lower triangle]
edgeWeights = [("l1", "l2", 1),("l1", "l2", 3),("l1", "l2", 4),("l1", "l2", 2),("l1", "l2", 2),("l1", "l2", 3),("l1", "l2", 3)]
mlvmst = ConstructMinLeavesVMST(edgeWeights)