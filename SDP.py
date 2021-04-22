#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 11:31:53 2021

@author: danielstilckfranca
"""

import time
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

import cvxgraphalgs as cvxgr
from cvxgraphalgs.structures import Cut

problem_graph = nx.random_regular_graph(d=4, n=12)

sdp_cut = cvxgr.algorithms.goemans_williamson_weighted(problem_graph)


nx.adjacency_matrix(problem_graph)