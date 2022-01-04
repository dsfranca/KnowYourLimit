
import numpy as np
import warnings
import cirq

from numpy import linalg as LA
from scipy import linalg


"""
The module containing the class we use to represent quantum channels.
"""

class quantum_channel:
    """
    The class for quantum channels. The quantum channel is initialized by a cirq noise model. We then compute the spectral gap and fixed point of the noise model.
    """
    def __init__(self, cirq_model):

        kraus_array=cirq.kraus(cirq_model)
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
        """
        Applies the quantum channel to a matrix.

        Returns:
            output (complex array): the output of the channel given an input matrix.

        Args:
    
            X(complex array): input of the channel.

        """
        output=np.zeros([2**self.locality,2**self.locality],np.complex128)
        for kraus_matrix in self.kraus:

            
            new=X.dot(kraus_matrix.getH())
            new=kraus_matrix.dot(new)
            output+=new
        return output
    
    #finds the overlap tr(A^*T(B)) for the channel T
    def overlap(self,A,B):
        """
        Denoting the underlying chanel by :math:`T` and given matrices :math:`A,B`, this function computes :math:`\\operatorname{tr}(A^{\\dagger}T(B))`.

        Returns:
            output (complex array): the output of the channel given an input matrix.

        Args:
    
            X(complex array): input of the channel.

        """
        new=self.apply_channel(B)
        new=(A.getH()).dot(new)
        return np.trace(new)
        





class lindbladian:
    """
    The class for Lindbladians. The quantum channel is initialized either by a cirq noise model or a process matrix in the Heisenberg picture. 
    We then compute the spectral gap and fixed point of the Lindbladian. It will srote the contraction rate towards the fixed point of the noise in self.contraction. The density matrix of the noise will be stored in self.fixed_point.



    Args:
    
        input_noise(complex array or cirq noise model): specification of the Lindbladian. This can either be in the form of a process matrix or a cirq noise model. In the latter case, the Lindbladian is assumed to be of the form :math:`\\mathcal{L}=r(T-\\textrm{id})`, where :math:`T` is the noise model and :math:`r` is the noise rate.

        rate (float): the rate :math:`r` of the noise as specified above. Only required if we are given a cirq noise model.

        cirq_input (bool): set to 1 if the input model for the noise model is a cirq quantum channel. If it is set to 0, it is assumed that we are given the process matrix of the Lindbladian.

    
    """
    def __init__(self, input_noise,rate=1,cirq_input=1):
        """
        The Lindbladian is initiated by a quantum channel.
        """
        if cirq_input is 1:
            qchannel=quantum_channel(input_noise)
            self.contraction=rate*qchannel.contraction
            self.fixed_point=qchannel.fixed_point
        else: 
            #compute spectrum
            self.transfer=input_noise
            [self.spectrum,self.eigenvectors]=linalg.eig(self.transfer)
            #compute singular values and spectrum
            self.U, self.s, self.Vh = linalg.svd(self.transfer)
            kernel_lind=linalg.null_space(self.transfer)
            
            #convert fixed point to quantum state
            self.fixed_point=np.reshape(kernel_lind[:,0],[2,2])
            self.fixed_point=self.fixed_point/np.trace(self.fixed_point)
            
            #w will now check that the conditions for the method to work apply to the noise model
            #check if classical
            check=self.fixed_point - np.diag(np.diag(self.fixed_point))
            if linalg.norm(check)>10**(-5):
                warnings.warn("Channel does not converge to classical state, method may yield overly optimistic estimates")

            if LA.matrix_rank(self.fixed_point, tol=10**-4)<self.dimension:
                warnings.warn("Fixed point is not of full-rank. Method does not apply.Consider approximating the noise by one with a fullrank.")
                
            
            if len(kernel_lind)>1:
                warnings.warn("Channel does not have unique fixed point. Method does not apply.")
            
            self.spectral_gap=np.sort(np.abs(self.spectrum))[-2]
            self.contraction=self.spectral_gap**2
        

            


