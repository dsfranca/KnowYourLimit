#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 10:43:35 2021

@author: danielstilckfranca
"""
from circuit_compiler import *
from partition_function_estimator import *
from noise_model import *
from circuit_analyser import *
import cirq
import warnings
import cvxgraphalgs as cvxgr
from cvxgraphalgs.structures import Cut
from scipy import sparse





def SDP_average(problem_graph,samples):
    results=[]
    for k in range(0,samples):
        sdp_cut = cvxgr.algorithms.goemans_williamson_weighted(problem_graph)
        A=nx.adjacency_matrix(problem_graph)
        results.append((0.5*np.sum(A)-2*sdp_cut.evaluate_cut_size(problem_graph)))
    return [np.mean(results),np.min(results)]
    
