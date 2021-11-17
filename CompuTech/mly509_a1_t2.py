"""
Based on odeint2.py, and attempts to solve task 2 of the assignment
"""

import matplotlib.pyplot as plt  # Plotting library

from scipy.integrate import odeint # SciPy ODE integration
from numpy import linspace, sqrt, exp

def task2(f, x, vs):
    # f is an array of all evolving variables
    phiHat = f[0]
    EHat = f[1]
    
    # Calculate the time derivatives
    dphiHat_by_dx = -EHat
    dEHat_by_dx = (1.0/sqrt(1.0-(2*phiHat/(vs*vs)))) - exp(phiHat)
    
    # Return the time derivatives in the same order as in the input f
    return [dphiHat_by_dx, dEHat_by_dx]


def run(vs):
    
    y0 = [0.0, 0.001] # Starting location 
    x = linspace(0, 40, 100)
    y = odeint(task2, y0, x, args = (vs,))

    plt.plot(x, y[:,0])
    plt.plot(x, y[:,1])
    plt.grid(True)  # Add a background grid
    plt.xlabel(r'$\hat{x}$ [normalised to Debye Length $λ_D$]')
    plt.title(r'Task2; $V_s$ = %.1f' % (vs))
    plt.legend([r'$\hat{\Phi} (\hat{x})$',r'$\hat{E}(\hat{x})$'])
    plt.annotate(r'$\hat{\Phi}$(40) = %.2f' % (y[100-1,0]),(40,y[100-1,0]), horizontalalignment='right')
    plt.annotate(r'$\hat{E}$(40) = %.2f' % (y[100-1,1]),(40,y[100-1,1]), horizontalalignment='right')
    plt.show()


if __name__ == "__main__":
    run(1.0)
    
