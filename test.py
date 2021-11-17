from scipy.integrate import odeint
from numpy import linspace, exp
import matplotlib.pyplot as plt
#Function to return dy/dt
def gradientFunc(curValue,curTime):
    return -10.*curValue
#Time for outputs
time=linspace(0,1,40)
y0= 10. #Initial condition
y=odeint(gradientFunc,y0,time) #Integrate
#Plot
plt.plot(time,y,'x',label='odeint')
plt.plot(time,y0*exp(-10.*time),label='Analytic')
plt.legend() ; plt.show()