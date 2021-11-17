'''
Assignment 1 Task 3, which depends on Task2 script
'''
from mly509_a1_t2 import task2
import matplotlib.pyplot as plt  # Plotting library
from scipy.integrate import odeint # SciPy ODE integration
from numpy import linspace, pi, sqrt, exp

if __name__ == "__main__":
    vs = 1
    y0 = [0.0, 0.001] # Starting location 
    x = linspace(0, 40, 100)
    y = odeint(task2, y0, x, args = (vs,))
    phiHat = y[:,0]
    jHat = sqrt(920.0/pi)*exp(phiHat)-1

    plt.plot(x, jHat)
    plt.grid(True)  # Add a background grid
    plt.xlabel(r'Distance from the plasma into the sheath [$\hat{x}$, normalised to Debye Length $λ_D$]')
    plt.ylabel(r'Current Density [Normalised to ion current density at wall $\hat{j}$]')
    plt.title(r'Task3; $V_s$ = %.1f' % (vs))
    plt.legend([r'$\hat{j} (\hat{x})$'])
    plt.ylim(-5,20)
    plt.xlim(0,40)
    plt.show()