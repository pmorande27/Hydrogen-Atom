from matplotlib import pyplot as plt
import numpy as np
import time
import scipy.interpolate
from scipy.sparse import linalg as splinalg
from scipy import linalg
from scipy.sparse import diags
from scipy import integrate
from scipy import optimize
from Hydrogen_Atom import Hydrogen_Atom
from Modified_Hydrogen_Atom import Hydrogen_Atom_M
c1 = 0.0380998 # nm^2 eV
c2 = 1.43996 # nm eV
r0 = 0.0529177 # nm
h  = 6.62606896e-34 # J s
c  = 299792458. # m/s
hc = 1239.8419 # eV nm
rmax = 10 # nm
def plot_error():
    """
    Function used to create a plot of the Error in the energy as a function of N,
    """
    N = [5*i for i in range(40,150)]
    Hs = [Hydrogen_Atom(n,rmax) for n in N]
    error = [H.get_error()[0] for H in Hs]
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
    H = Hydrogen_Atom(N,rmax)
    er1 ,er2= H.get_error()
    while er1> 5e-4 or  er2 > 5e-4:
        N = N+1
        er1 ,er2= H.get_error()
    return N
def energy_dif_vs_alpha_plot(N=1024):
    """
    This funcntion makes the Energy dfference vs alpha plot requested in the exercise for a given N, the defaul is assumed
    to be 1024. The range of values of alpha is not specified so it is taken to be (0,0.99) inclusive.
    """
    alpha_linspace = np.linspace(0,0.99,50)
    Hms = [Hydrogen_Atom_M(N,alpha) for alpha in alpha_linspace]
    diffs = [Hm.get_energy_difference_super_modified() for Hm in Hms]
    alpha = [0,0.01]
    plt.plot(alpha_linspace,diffs)
    plt.xlabel('α')
    plt.ylabel('ΔE[eV]')
    plt.show()
#energy_dif_vs_alpha_plot()  
def plot_my_v(N=1024):
    """
    Function used to create a plot of the new modified potential with alpha = 0.01.
    N will give the number of points and it will explore the values between rmax/N and rmax (by choice as it is not specified)
    as they are relevant for the following sections.
    """
    r = np.linspace(rmax/N, rmax, N)
    alpha = 0.01
    Hm = Hydrogen_Atom_M(N,alpha)
    v_my = [Hm.potential_numerical(x) for x in r]
    plt.plot(r,v_my)
    plt.xlabel('r[nm]')
    plt.ylabel('V[eV]')
    plt.show()  
#plot_my_v(100)

def plot_lambda_error():
    """
    Function used to get a plot of how the error in lambda (the difference between calculated and expected values) evolves
    with different alphas
    """
    N = 1024
    alphas = np.linspace(0,0.01,100)
    lambda_r = 121.5
    Hms = [Hydrogen_Atom_M(N,alpha) for alpha in alphas]
    error = [abs(hc/Hm.get_energy_difference_super_modified() - lambda_r) for Hm  in Hms]
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
    Hm = Hydrogen_Atom_M(N,alpha)
    delta_e = Hm.get_energy_difference_super_modified()
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
        Hm = Hydrogen_Atom_M(N,alpha)
        delta_e = Hm.get_energy_difference_super_modified()
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

def main():
    '''
    N = 636
    H = Hydrogen_Atom(N)
    er1 ,er2= H.get_error()
    print (f"Err1 = {er1}, Err2 = {er2}.")
    assert er1 < 5e-4
    assert er1 < 5e-4
    print(H.get_energy_levels()) 
    N = 1204
    Hm = Hydrogen_Atom_M(N,alpha=0.01)
    E1, E2 = Hm.energy_levels_modified_super()
    E1_th, E2_th = -13.807387841665346, -3.5346025272551795
    Err1 = abs((E1 - E1_th) / E1_th)
    Err2 = abs((E2 - E2_th) / E2_th)
    print (f"{Err1 = }, {Err2 = }")
    assert Err1 < 5e-4
    assert Err2 < 5e-4
    '''
    r = 1000 # nm
    N = 100000
    number = 10
    H = Hydrogen_Atom(N,r)
    t1 = time.time()
    vals = H.calculate_energy_levels_super(number)
    my_e1 = vals[0]
    my_e2 = vals[1]
    #print(vals)
    t2 = time.time()
    print (f"Calculation took {t2-t1} seconds.")
    theo = [-c2/(2*r0*n**2) for n in range(1,101)]
    error = [abs((x-y)/x) for (x,y) in zip(theo,vals)]
    print(error)
    e1_th = -c2 / (2 * r0)
    e2_th = e1_th / 4

    er1 = abs((my_e1 - e1_th) / e1_th)
    er2 = abs((my_e2 - e2_th) / e2_th)
    print (f"Err1 = {er1}, Err2 = {er2}.")
    plt.plot(vals,'o')
    plt.show()

main()