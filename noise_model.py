#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 18:18:21 2021

@author: danielstilckfranca
"""
import numpy as np
import warnings
import cirq

from numpy import linalg as LA
from scipy import linalg


class quantum_channel:
    def __init__(self, cirq_model):
        
        kraus_array=cirq_model
        #converts representation from ndarray to matrix, easier for complex conjugate
        self.kraus=[]
        for kraus in kraus_array:
            self.kraus.append(np.matrix(kraus))
            
        tuple(self.kraus)
        self.locality=int(np.log(np.shape(self.kraus[0])[0])/np.log(2)) 
        self.dimension=2**self.locality
        #find the matrix representation of the channel in the standard matrix basis
        self.transfer=np.matrix(np.zeros([self.dimension**2,self.dimension**2],dtype=np.complex128))
        for i in range(0,self.dimension):
            for j in range(0,self.dimension):
                for k in range(0,self.dimension):
                    for l in range(0,self.dimension):
                        left_vector=np.matrix(np.zeros([self.dimension,self.dimension],dtype=np.complex128))
                        left_vector[i,j]=1
                        right_vector=np.matrix(np.zeros([self.dimension,self.dimension],dtype=np.complex128))
                        right_vector[k,l]=1
                        
                        self.transfer[i*self.dimension+j,k*self.dimension+l]=self.overlap(left_vector,right_vector)
        
        #compute spectrum
        [self.spectrum,self.eigenvectors]=linalg.eig(self.transfer)
        #compute singular values and spectrum
        self.U, self.s, self.Vh = linalg.svd(self.transfer)
        
        
        #convert fixed point to quantum state
        self.fixed_point=np.reshape(self.eigenvectors[:,0],[self.dimension,self.dimension])
        self.fixed_point=self.fixed_point/np.trace(self.fixed_point)
        
        #w will now check that the conditions for the method to work apply to the noise model
        #check if classical
        check=self.fixed_point - np.diag(np.diag(self.fixed_point))
        if linalg.norm(check)>10**(-5):
            warnings.warn("Channel does not converge to classical state, method may yield overly optimistic estimates")

        if LA.matrix_rank(self.fixed_point, tol=10**-4)<self.dimension:
            warnings.warn("Fixed point is not of full-rank. Method does not apply.Consider approximating the noise by one with a fullrank.")
            
        
        if np.abs(self.spectrum[1])>=1-10**(-7):
            warnings.warn("Channel does not have unique fixed point. Method does not apply.")
        
        self.spectral_gap=np.sort(np.abs(self.spectrum))[1]
        self.contraction=self.spectral_gap**2
    #applies the channel to a matrix X    
    def apply_channel(self,X):
        output=np.zeros([2**self.locality,2**self.locality],np.complex128)
        for kraus_matrix in self.kraus:

            
            new=X.dot(kraus_matrix.getH())
            new=kraus_matrix.dot(new)
            output+=new
        return output
    
    #finds the overlap tr(A^*T(B)) for the channel T
    def overlap(self,A,B):
        new=self.apply_channel(B)
        new=(A.getH()).dot(new)
        return np.trace(new)
        








