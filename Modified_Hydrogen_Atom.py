from matplotlib import pyplot as plt
import numpy as np
import time
import scipy.interpolate
from scipy.sparse import linalg as splinalg
from scipy import linalg
from scipy.sparse import diags
from scipy import integrate
from scipy import optimize
c1 = 0.0380998 # nm^2 eV
c2 = 1.43996 # nm eV
r0 = 0.0529177 # nm
h  = 6.62606896e-34 # J s
c  = 299792458. # m/s
hc = 1239.8419 # eV nm
rmax = 1.5 # nm
class Hydrogen_Atom_M(object):
    def __init__(self,N,alpha) -> None:
        self.N = N
        self.alpha = alpha

    def potential_numerical(self,r):
        """
        Function used to calculate the integral of the potential Energy
        """
        y = lambda x: -np.power(x,self.alpha-2)*c2/np.power(r0,self.alpha)
        return integrate.quad(y,r,np.inf)[0]
    def potential_exact(self,r):
        return c2*np.power(r,self.alpha-1)*np.power(r0,-self.alpha) / (self.alpha-1)

    def laplacian_operator(self):
        """
        This function creates an N by N matrix that represents the laplacian operator in one dimension,
        the range of values of x explored will be spaced by rmax/N (tridiagonal matrix), it will return the laplacian
        multiplied by -c1 so that the construction of the hamiltonian is easier
        """
        h = rmax/self.N
        laplacian = diags([-c1/np.power(h,2), 2*c1/np.power(h,2), -c1/np.power(h,2)], [-1, 0, 1], shape=(self.N, self.N))
        return(laplacian)



    def potential_operator_modified(self):
        """
        Used to get the matrix representation (N x N) of the new modified operator, it depends on alpha
        """
        r_space = np.linspace(rmax/self.N,rmax,self.N)
        diagonals =  np.zeros(self.N)
        for i in range(self.N):
            diagonals[i] = self.potential_numerical(r_space[i],self.alpha)
        diagonals = [diagonals]
        potential_matrix = diags(diagonals, [0])
        return potential_matrix
    def calculate_energy_levels(self):
        """
        Used to calculate the new energy levels (the first two eigenvalues)
        """
        hamitltonian_operator= self.laplacian_operator()+self.potential_operator_modified()
        vals, vecs = splinalg.eigsh(hamitltonian_operator, k=2,which= 'SA')
        return vals
    '''
    N = 1024
    E1, E2 = calculate_energy_levels(N, alpha=0.01)
    E1_th, E2_th = -13.807387841665346, -3.5346025272551795
    Err1 = abs((E1 - E1_th) / E1_th)
    Err2 = abs((E2 - E2_th) / E2_th)
    print (f"{Err1 = }, {Err2 = }")
    assert Err1 < 5e-4
    assert Err2 < 5e-4
    '''
    def get_energy_difference(self):
        """
        Gets the difference in energy of the two obtained levels using the new hamiltonian
        """
        E1, E2 = self.calculate_energy_levels()
        return (E2-E1)
    def energy_levels_modified_super(self):
        """
        Used to calculate the new energy levels (the first two eigenvalues) but faster using the result of task 5
        it will use the advantage that H is tridiagonal and Hermitian to speed things up
        """
        h = np.divide(rmax,self.N)
        laplacian = diags([-c1/np.power(h,2), 2*c1/np.power(h,2), -c1/np.power(h,2)], [-1, 0, 1], shape=(self.N, self.N))
        diagonal_off = -c1/np.power(h,2)*np.ones(self.N-1)
        r_space = np.linspace(rmax/self.N,rmax,self.N)
        diagonal_potential =  np.zeros(self.N)
        for i in range(self.N):
            diagonal_potential[i] =self.potential_numerical(r_space[i])
        diagonal_main =  2*c1/np.power(h,2)*np.ones(self.N)+diagonal_potential
        #print(diagonal_off)
        #print(diagonal_main)
        val_one, val_two = linalg.eigvalsh_tridiagonal(diagonal_main,diagonal_off,select='i',select_range=(0,1))
        return val_one,val_two
    '''
    N = 1204
    E1, E2 = energy_levels_modified_super(N, alpha=0.01)
    E1_th, E2_th = -13.807387841665346, -3.5346025272551795
    Err1 = abs((E1 - E1_th) / E1_th)
    Err2 = abs((E2 - E2_th) / E2_th)
    print (f"{Err1 = }, {Err2 = }")
    assert Err1 < 5e-4
    assert Err2 < 5e-4
    '''
    def get_energy_difference_super_modified(self):
        """
        Gets the difference in energy of the two obtained levels using the new hamiltonian and using the fastest method
        """
        E1, E2 = self.energy_levels_modified_super()
        return (E2-E1)
