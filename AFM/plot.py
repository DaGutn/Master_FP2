import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pylab import *

Spurabstand = ufloat (1.52,0.002)*10**(-6)
Pitlaenge = ufloat (0.78,0.02)*10**(-6)

print ((4/(Spurabstand*Pitlaenge))*0.0089/17)