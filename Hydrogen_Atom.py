## Imports 
from matplotlib import pyplot as plt
import numpy as np
import time
import scipy.interpolate
from scipy.sparse import linalg as splinalg
from scipy import linalg
from scipy.sparse import diags
from scipy import integrate
from scipy import optimize
# Constants (use these)
# Constants (use these)
c1 = 0.0380998 # nm^2 eV
c2 = 1.43996 # nm eV
r0 = 0.0529177 # nm
h  = 6.62606896e-34 # J s
c  = 299792458. # m/s
hc = 1239.8419 # eV nm
class Hydrogen_Atom(object):
    def __init__(self,N,rmax) -> None:
        self.rmax = rmax
        self.N = N
    def laplacian_operator(self):
        """
        This function creates an N by N matrix that represents the laplacian operator in one dimension,
        the range of values of x explored will be spaced by rmax/N (tridiagonal matrix), it will return the laplacian
        multiplied by -c1 so that the construction of the hamiltonian is easier
        """
        h = self.rmax/self.N
        laplacian = diags([-c1/np.power(h,2), 2*c1/np.power(h,2), -c1/np.power(h,2)], [-1, 0, 1], shape=(self.N, self.N))
        return(laplacian)
    def potential_operator(self):
        """
        Function used to create an N by N matrix that represent the potential operator,
        it will be a diagonal matrix.
        """
        
        r_space = np.linspace(self.rmax/self.N,self.rmax,self.N)
        diagonals = [-c2/r_space]
        potential_matrix = diags(diagonals, [0])
        
        return potential_matrix

    def hamiltonian_operator(self):
        """
        Function used to create an N by N matrix that represent the Hamiltonian operator,
        it will be a diagonal matrix.
        """
        return self.laplacian_operator()+self.potential_operator()
    def get_energy_levels(self):
        """
        Function used to get the first two eigenvalues of the hamiltonian which correspond to the Energy levels of the system.
        """
        vals, vecs = splinalg.eigsh(self.hamiltonian_operator(), k=10,which= 'SA')
        print(vals)
        return vals[0].real,vals[1].real

    def get_error(self):
        """
        Function used to get the error in the energy levels when considering an N by N hamiltonian
        """
        e_calc_one,e_calc_two = self.get_energy_levels()
        e_theoretical_one = -c2*1/(2*r0)
        e_theoretical_two = -c2*1/(8*r0)
        e_error_one = abs((e_calc_one-e_theoretical_one)/e_theoretical_one)
        e_error_two = abs((e_calc_two-e_theoretical_two)/e_theoretical_two)
        return e_error_one,e_error_two
    def calculate_energy_levels_super(self):
        """
        Function used to calculate the energy levels but in a faster way, it uses the fact that H is both tridiagonal and 
        hermitian
        """
        t1 = time.time()
        h = np.divide(self.rmax,self.N)
        laplacian = diags([-c1/np.power(h,2), 2*c1/np.power(h,2), -c1/np.power(h,2)], [-1, 0, 1], shape=(self.N, self.N))
        diagonal_off = -c1/np.power(h,2)*np.ones(self.N-1)
        r_space = np.linspace(self.rmax/self.N,self.rmax,self.N)
        diagonal_potential = -c2/r_space
        diagonal_main =  2*c1/np.power(h,2)*np.ones(self.N)+diagonal_potential
        val_one, val_two = linalg.eigvalsh_tridiagonal(diagonal_main,diagonal_off,select='i',select_range=(0,1))
        vals  = linalg.eigvalsh_tridiagonal(diagonal_main,diagonal_off,select='i',select_range=(0,10))
        print(vals)

        return val_one,val_two

'''
N = 636
er1 ,er2= get_error(N)
print (f"Err1 = {er1}, Err2 = {er2}.")
assert er1 < 5e-4
assert er1 < 5e-4
print(get_energy_levels(636))
'''
'''
t1 = time.time()
my_e1, my_e2 = calculate_energy_levels_super()
t2 = time.time()
print (f"Calculation took {t2-t1} seconds.")

e1_th = -c2 / (2 * r0)
e2_th = e1_th / 4

er1 = abs((my_e1 - e1_th) / e1_th)
er2 = abs((my_e2 - e2_th) / e2_th)
print (f"Err1 = {er1}, Err2 = {er2}.")
'''
