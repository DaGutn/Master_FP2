import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pylab import *


x0_1=ufloat (641.771,0.018)*10**(-9)
x0_2=ufloat (581.632,0.014)*10**(-9)
x0_3=ufloat (544.061,0.015)*10**(-9)
x0_4=ufloat (539.09,0.14)*10**(-9)
me=0.13*9.11*10**(-31)
mh=0.45*9.11*10**(-31)

mu=0.13*(1.602*10**(-19))
#mu=(me*mh)/(me+mh)
pi=3.141592653589793

Er=9.5
Eg=1.74*(1.602*10**(-19))
e=1.602*10**(-19)
E0=8.85*10**(-12)
h=6.626*10**(-34)/(2*pi)

A=((pi**2*h*h)/2)*(1/me+1/mh)
B=(mu*e**4)/(32*pi**2*E0**2*Er**2*h**2)
C=(1.786*e**2)/(4*pi*Er*E0)
D=(h**2*pi**2)/(2*mu)

#x0_1=(6.626*10**(-34)*3*10**8)/x0_1
print (x0_1)

#a1=1*C/(2*(Eg-x0_1-B))+((C/(2*(Eg-x0_1-B)))**2-(A/(Eg-x0_1-B)))**(1/2)
#a2=1*C/(2*(Eg-x0_2-B))+((C/(2*(Eg-x0_2-B)))**2-(A/(Eg-x0_2-B)))**(1/2)

a1=C/(2*(Eg-x0_1))+((C/(2*(Eg-x0_1)))**2-(D/(Eg-x0_1)))**(1/2)

print(a1)
#print(a2)

#Polarisation
Pol0_1=ufloat(0.61438,0.000728)
Pol0_2=ufloat(0.7097,0.000758)
Pol0_3=ufloat(0.86823,0.000695)
Pol0_4=ufloat(1.0043,0.00297)

Pol90_1=ufloat(0.61438,0.000144)
Pol90_2=ufloat(0.70561,0.00075)
Pol90_3=ufloat(0.86621,0.000696)
Pol90_4=ufloat(0.99321,0.00295)

P_1=(Pol0_1-Pol90_1)/(Pol0_1+Pol90_1)
P_2=(Pol0_2-Pol90_2)/(Pol0_2+Pol90_2)
P_3=(Pol0_3-Pol90_3)/(Pol0_3+Pol90_3)
P_4=(Pol0_4-Pol90_4)/(Pol0_4+Pol90_4)

#print(P_1)
#print(P_2)
#print(P_3)
#print(P_4)




Spurabstand = ufloat (1.52,0.002)*10**(-6)
Pitlaenge = ufloat (0.78,0.02)*10**(-6)

#print ((4/(Spurabstand*Pitlaenge))*0.0089/17)