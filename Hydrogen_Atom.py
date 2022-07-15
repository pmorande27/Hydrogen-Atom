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
rmax = 1.5 # nm
def laplacian_operator(N):
    """
    This function creates an N by N matrix that represents the laplacian operator in one dimension,
    the range of values of x explored will be spaced by rmax/N (tridiagonal matrix), it will return the laplacian
    multiplied by -c1 so that the construction of the hamiltonian is easier
    """
    h = rmax/N
    laplacian = diags([-c1/np.power(h,2), 2*c1/np.power(h,2), -c1/np.power(h,2)], [-1, 0, 1], shape=(N, N))
    return(laplacian)
def potential_operator(N):
    """
    Function used to create an N by N matrix that represent the potential operator,
    it will be a diagonal matrix.
    """
    
    r_space = np.linspace(rmax/N,rmax,N)
    diagonals = [-c2/r_space]
    potential_matrix = diags(diagonals, [0])
    
    return potential_matrix

def hamiltonian_operator(N):
    """
    Function used to create an N by N matrix that represent the Hamiltonian operator,
    it will be a diagonal matrix.
    """
    return laplacian_operator(N)+potential_operator(N)
def get_energy_levels(N):
    """
    Function used to get the first two eigenvalues of the hamiltonian which correspond to the Energy levels of the system.
    """
    vals, vecs = splinalg.eigsh(hamiltonian_operator(N), k=2,which= 'SA')
    return vals[0].real,vals[1].real

def get_error(N):
    """
    Function used to get the error in the energy levels when considering an N by N hamiltonian
    """
    e_calc_one,e_calc_two = get_energy_levels(N)
    e_theoretical_one = -c2*1/(2*r0)
    e_theoretical_two = -c2*1/(8*r0)
    e_error_one = abs((e_calc_one-e_theoretical_one)/e_theoretical_one)
    e_error_two = abs((e_calc_two-e_theoretical_two)/e_theoretical_two)
    return e_error_one,e_error_two
def plot_error():
    """
    Function used to create a plot of the Error in the energy as a function of N,
    """
    N = [5*i for i in range(40,150)]
    error = [get_error(n)[0] for n in N]
    plt.plot(N,error)
    plt.xlabel("N")
    plt.ylabel("error in energy [eV]")
    plt.show()
#plot_error()
# To generate the graph un-comment the line above.
def get_minimum_N():
    """
    Function used to get the minimum value of N that gives a precission of 5e-4 (0.05%) to both energy levels
    """
    N = 600
    er1 ,er2= get_error(N)
    while er1> 5e-4 or  er2 > 5e-4:
        N = N+1
        er1 ,er2= get_error(N)
    return N
'''
N = 636
er1 ,er2= get_error(N)
print (f"Err1 = {er1}, Err2 = {er2}.")
assert er1 < 5e-4
assert er1 < 5e-4
print(get_energy_levels(636))
'''
def potential_numerical(r, alpha):
    """
    Function used to calculate the integral of the potential Energy
    """
    y = lambda x: -np.power(x,alpha-2)*c2/np.power(r0,alpha)
    return integrate.quad(y,r,np.inf)[0]
def potential_exact(r, alpha):
    return c2*np.power(r,alpha-1)*np.power(r0,-alpha) / (alpha-1)
'''
for my_r in np.linspace(0.01, 1, 100):
    diff = abs(potential_numerical(my_r, 0.01) - potential_exact(my_r, 0.01))
    assert(diff <= 1e-5)
'''


def plot_my_v(N=1024):
    """
    Function used to create a plot of the new modified potential with alpha = 0.01.
    N will give the number of points and it will explore the values between rmax/N and rmax (by choice as it is not specified)
    as they are relevant for the following sections.
    """
    r = np.linspace(rmax/N, rmax, N)
    alpha = 0.01
    v_my = [potential_numerical(x,alpha)for x in r]
    plt.plot(r,v_my)
    plt.xlabel('r[nm]')
    plt.ylabel('V[eV]')
    plt.show()
plot_my_v(100)

alpha = 0.01
def potential_operator_modified(N,alpha):
    """
    Used to get the matrix representation (N x N) of the new modified operator, it depends on alpha
    """
    r_space = np.linspace(rmax/N,rmax,N)
    diagonals =  np.zeros(N)
    for i in range(N):
        diagonals[i] = potential_numerical(r_space[i],alpha)
    diagonals = [diagonals]
    potential_matrix = diags(diagonals, [0])
    return potential_matrix
def calculate_energy_levels(N,alpha):
    """
    Used to calculate the new energy levels (the first two eigenvalues)
    """
    hamitltonian_operator= laplacian_operator(N)+potential_operator_modified(N,alpha)
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
def get_energy_difference(N,alpha):
    """
    Gets the difference in energy of the two obtained levels using the new hamiltonian
    """
    E1, E2 = calculate_energy_levels(N, alpha)
    return (E2-E1)
def energy_levels_modified_super(N,alpha):
    """
    Used to calculate the new energy levels (the first two eigenvalues) but faster using the result of task 5
    it will use the advantage that H is tridiagonal and Hermitian to speed things up
    """
    h = np.divide(rmax,N)
    laplacian = diags([-c1/np.power(h,2), 2*c1/np.power(h,2), -c1/np.power(h,2)], [-1, 0, 1], shape=(N, N))
    diagonal_off = -c1/np.power(h,2)*np.ones(N-1)
    r_space = np.linspace(rmax/N,rmax,N)
    diagonal_potential =  np.zeros(N)
    for i in range(N):
        diagonal_potential[i] = potential_numerical(r_space[i],alpha)
    diagonal_main =  2*c1/np.power(h,2)*np.ones(N)+diagonal_potential
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
def get_energy_difference_super_modified(N,alpha):
    """
    Gets the difference in energy of the two obtained levels using the new hamiltonian and using the fastest method
    """
    E1, E2 = energy_levels_modified_super(N, alpha)
    return (E2-E1)
def energy_dif_vs_alpha_plot(N=1024):
    """
    This funcntion makes the Energy dfference vs alpha plot requested in the exercise for a given N, the defaul is assumed
    to be 1024. The range of values of alpha is not specified so it is taken to be (0,0.99) inclusive.
    """
    alpha_linspace = np.linspace(0,0.99,50)
    diffs = [get_energy_difference_super_modified(N,x) for x in alpha_linspace]
    alpha = [0,0.01]
    plt.plot(alpha_linspace,diffs)
    plt.xlabel('α')
    plt.ylabel('ΔE[eV]')
    plt.show()
#energy_dif_vs_alpha_plot()    
def plot_lambda_error():
    """
    Function used to get a plot of how the error in lambda (the difference between calculated and expected values) evolves
    with different alphas
    """
    N = 1024
    alphas = np.linspace(0,0.01,100)
    lambda_r = 121.5
    error = [abs(hc/get_energy_difference_super_modified(N,alpha) - lambda_r) for alpha in alphas]
    plt.plot(alphas,error)
    plt.xlabel("α")
    plt.ylabel("error[nm]")
    plt.show()
#plot_lambda_error()
# To generate the plot un-comment the line above, it can take a while
def get_lambda(alpha,epsilon,N):
    """
    Function used to represent the equation |lambda_obtained - lambda_real| -0.1 + epsilon, the root of this equation will
    give a close approximation of what is alpha max, epsilon is a correction parameter initially set to 0, if the given result
    of alpha results in an error bigger than 0.1 epsilon will increase by a fixed and small amount so that the algorithm finds
    another root close to the first one but more leaning towards an error smaller than 0.1
    """
    lambda_r = 121.5
    delta_e = get_energy_difference_super_modified(N,alpha)
    lambda_ = hc/delta_e
    return abs(lambda_-lambda_r)-0.1+epsilon
def find_alpha_max(N = 1024):
    """
    Function used to get the maximum alpha that produces an error less than 0.1, for that it starts by calculating an 
    approximation of the root for the equation |lambda_obtained - lambda_real| -0.1 = 0, if the 
    calculated value of alpha produces an error smaller than 0.1 then the result is accepted,
    but if the result calculated gives an error bigger than 0.1 then it tries to use the epsilo parameter to find a solution
    of alpha that is close to that initial guess but that produces a smaller error.
    """
    epsilon = 0
    step = 10e-8    
    alpha = scipy.optimize.brentq(lambda x: get_lambda(x,epsilon,N),0.001,0.002,xtol=2e-13)
    while True:
        delta_e = get_energy_difference_super_modified(N,alpha)
        lambda_r = 121.5
        lambda_ = hc/delta_e
        error = abs(lambda_-lambda_r)
        if error < 0.1:
            break
        else:
            epsilon +=step
            alpha = scipy.optimize.brentq(lambda x: get_lambda(x,epsilon,N),0.001,0.002,xtol=2e-13)
    print(error)
    return alpha

'''
amax = find_alpha_max()
print (f"alpha_max = {amax}.")
'''
def plot_N_dependence():
    """
    Function used to show the dependence of alpha_max of the value of N
    """
    exponents = [1024+100]
    N = [1024 +200*i for i in range(30)]
    amaxs = [find_alpha_max(n) for n in N]
    plt.plot(N,amaxs)
    plt.xlabel('N')
    plt.ylabel('α_max')
    plt.show()
#plot_N_dependence()
# To generate the plot un-comment the line above, it can take a while but the plot should be obvious enough
# To get an approximation of a better value of alpha run this cell, it should take about 2 minutes
#amax = find_alpha_max(9000)
#print (f"alpha_max = {amax}.")
def calculate_energy_levels_super(N=100000):
    """
    Function used to calculate the energy levels but in a faster way, it uses the fact that H is both tridiagonal and 
    hermitian
    """
    t1 = time.time()
    h = np.divide(rmax,N)
    laplacian = diags([-c1/np.power(h,2), 2*c1/np.power(h,2), -c1/np.power(h,2)], [-1, 0, 1], shape=(N, N))
    diagonal_off = -c1/np.power(h,2)*np.ones(N-1)
    r_space = np.linspace(rmax/N,rmax,N)
    diagonal_potential = -c2/r_space
    diagonal_main =  2*c1/np.power(h,2)*np.ones(N)+diagonal_potential
    val_one, val_two = linalg.eigvalsh_tridiagonal(diagonal_main,diagonal_off,select='i',select_range=(0,1))
    return val_one,val_two
t1 = time.time()
my_e1, my_e2 = calculate_energy_levels_super()
t2 = time.time()
print (f"Calculation took {t2-t1} seconds.")

e1_th = -c2 / (2 * r0)
e2_th = e1_th / 4

er1 = abs((my_e1 - e1_th) / e1_th)
er2 = abs((my_e2 - e2_th) / e2_th)
print (f"Err1 = {er1}, Err2 = {er2}.")