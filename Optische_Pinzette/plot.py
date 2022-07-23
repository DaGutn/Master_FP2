from cProfile import label
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import scipy.constants as const


def laserleistung(A):
    B = np.zeros(len(A))
    for i in np.arange(len(B)):
        if(38 < A[i] <= 60):
            B[i] = 0.029*A[i]**2 - 2.139*A[i] + 39.23
        elif(60 < A[i] <= 350):
            B[i] = 0.598*A[i] - 26.3
    return B

### Aufgabe 1
d_kugel = 2.06 * 10**(-6)
pixelgroesse = d_kugel/122

### Aufgabe 2
## b)
current, m_x, m_y = np.genfromtxt('data/Konversionsfaktoren.txt', skip_header=1, unpack=True)
power = laserleistung(current)

fig, ax = plt.subplots(figsize=(4.5,3))
ax.plot(power, m_x, 'x--', label='x-Konversionsfaktor')
ax.plot(power, m_y, 'x--', label='y-Konversionsfaktor')

ax.set_xlabel(r'$P$ [mW]')
ax.set_ylabel(r'$m_{\mathrm{konv}}$ [V/$\mu$m]')
ax.legend(loc='best')

plt.tight_layout()
plt.savefig('plots/Konversionsfaktoren.pdf')
# plt.show()
plt.close()

m_x_mean = np.abs(np.mean(m_x))
m_y_mean = np.abs(np.mean(m_y))

## c)
z, voltage = np.genfromtxt('data/Diodensumme_zPiezo.txt', skip_header=1, unpack=True)

# ####################################
# def I_z(z, z_0):
#     return np.sqrt(1 + (z/z_0)**2) * np.sin(np.arctan(z/z_0))

# lam = 975*10**(-9)
# n=1.33
# w_0 = 1.22*lam/n * np.sqrt((n/1.25)**2 - 1)
# z_0 = np.pi*w_0**2/lam *10**(7)

# params, cov = curve_fit(I_z, z[15:23], z_0)
# z_fit = np.linspace(28.5, 34)
# ####################################

fig, ax = plt.subplots(figsize=(4.5,3))
ax.plot(z, voltage, 'x--', color='green')
# ax.plot(z_fit, I_z(z_fit, *params))
ax.axvline(33, color='grey', linestyle='dashed', linewidth=2, alpha=0.8, label='Probenoberflächenposition')

ax.set_xlabel(r'$z$ [$\mu$m]')
ax.set_ylabel(r'Diodenspannung [V]')
ax.legend(loc='best')

plt.tight_layout()
plt.savefig('plots/Diodensumme.pdf')
# plt.show()
plt.close()

###################################################################################################
### Aufgabe 3 Fallensteifigkeit, Boltzmannkonstante
def PSD(f, A, f_0):
    return A/(f**2 + f_0**2)

def computeK(f_0):
    return 3 * np.pi * d_kugel * 8.9*10**(-4) * 2 * np.pi *f_0

def computeVar(x):
    return np.var(x - np.mean(x))

def computeK_B(k, x):
    return computeVar(x) * k / (273.15 + 25)

def linear(x, a, b):
    return a*x + b


# ## a) Fallensteifigkeit als Funktion der Laserleistung
# current = [150, 200, 250, 300, 350]
# power = laserleistung(current)
# print(power)
# frequenz = []
# PSD_x = []
# PSD_y = []

# for i in np.arange(len(current)):
#     with open('data/k_keineKraft_' + str(current[i]) + 'mA.FDdat', 'rb') as f:
#         clean_lines = (line.replace(b',',b'') for line in f)
#         frequenz_temp, PSD_x_temp, PSD_y_temp, egal = np.genfromtxt(clean_lines, unpack=True)
#         frequenz.append(frequenz_temp)
#         PSD_x.append(PSD_x_temp)
#         PSD_y.append(PSD_y_temp)

# f_min = [5, 5, 5, 5, 5]
# f_max = [50000, 50000, 50000, 50000, 50000]
# PSD_x_params = np.zeros((len(current), 2))
# PSD_y_params = np.zeros((len(current), 2))
# for i in np.arange(len(current)):
#     PSD_x_params[i], _ = curve_fit(PSD, frequenz[i][f_min[i]:f_max[i]], PSD_x[i][f_min[i]:f_max[i]])
#     PSD_y_params[i], _ = curve_fit(PSD, frequenz[i][f_min[i]:f_max[i]], PSD_y[i][f_min[i]:f_max[i]])

# # num = 3
# # frequenz_fit = np.linspace(0, frequenz[num][-1], 10000)
# # fig, axs = plt.subplots(2, 1, figsize=(6.5,5))

# # axs[0].plot(frequenz[num], PSD_x[num], label=r'Daten PSD$_x$')
# # axs[0].plot(frequenz_fit, PSD(frequenz_fit, *PSD_x_params[num]), label='Lorentzfunktion-Fit')
# # axs[0].set_xscale('log')
# # axs[0].set_yscale('log')
# # axs[0].set_xlim(right=6e4)
# # axs[0].set_xlabel(r'$f$ [Hz]')
# # axs[0].set_ylabel(r'PSD$_x$ [ms$^{\frac{1}{2}}$]')
# # axs[0].legend(loc='best')

# # axs[1].plot(frequenz[num], PSD_y[num], label=r'Daten PSD$_y$')
# # axs[1].plot(frequenz_fit, PSD(frequenz_fit, *PSD_y_params[num]), label='Lorentzfunktion-Fit')
# # axs[1].set_xscale('log')
# # axs[1].set_yscale('log')
# # axs[1].set_xlim(right=6e4)
# # axs[1].set_xlabel(r'$f$ [Hz]')
# # axs[1].set_ylabel(r'PSD$_y$ [ms$^{\frac{1}{2}}$]')
# # axs[1].legend(loc='best')

# # plt.tight_layout()
# # plt.savefig('plots/PSD_keineKraft.pdf')
# # # plt.show()
# # plt.close()


# # Fallensteifigkeit:
# k_x = computeK(PSD_x_params[:,1])
# k_y = computeK(PSD_y_params[:,1])

# print('================== Ohne Kraft:')
# for i in np.arange(len(current)):
#     print('x- und y-Roll-Off-Frequenz für ' + str(current[i]) + 'mA: {0:.3e}     {1:.3e}'.format(PSD_x_params[i,1], PSD_y_params[i,1]))

# for i in np.arange(len(current)):
#     print('x- und y-Fallensteifigkeit für ' + str(current[i]) + 'mA: {0:.3e}     {1:.3e}'.format(k_x[i], k_y[i]))

# k_x_params, _ = curve_fit(linear, power, k_x)
# k_y_params, _ = curve_fit(linear, power, k_y)

# # fig, ax = plt.subplots(figsize=(5.5,3.5))
# # ax.plot(power, k_x, 'x', label=r'Daten $k_x$')
# # ax.plot(power, k_y, 'x', label=r'Daten $k_y$')
# # power_fit = np.linspace(50, 200)
# # ax.plot(power_fit, linear(power_fit, *k_x_params), c='C0', ls='dashed', label='Linearer Fit')
# # ax.plot(power_fit, linear(power_fit, *k_y_params), c='C1', ls='dashed', label='Linearer Fit')

# # ax.set_xlabel(r'$P$ [mW]')
# # ax.set_ylabel(r'$k$ [N/m]')
# # ax.legend(loc='best')

# # plt.tight_layout()
# # plt.savefig('plots/k_keineKraft.pdf')
# # # plt.show()
# # plt.close()

# # Boltzmann-Konstante
# t = []
# x = []
# y = []
# for i in np.arange(len(current)):
#     t_temp, x_temp, y_temp, egal = np.genfromtxt('data/k_keineKraft_' + str(current[i]) + 'mA.TDdat', unpack=True)
#     t.append(t_temp)
#     x.append(x_temp)
#     y.append(y_temp)

# x = x / m_x_mean
# y = y / m_y_mean

# k_B_x = computeK_B(k_x, x*10**(-6))
# k_B_y = computeK_B(k_y, y*10**(-6))

# for i in np.arange(len(current)):
#     print('x- und y-Boltzmannkonst. für ' + str(current[i]) + 'mA: {0:.3e}     {1:.3e}'.format(k_B_x[i], k_B_y[i]))

# for i in np.arange(len(current)):
#     print('x- und y-Abweichungen für ' + str(current[i]) + 'mA: {0:.3f}     {1:.3f}'.format((k_B_x[i]-const.k)/const.k, (k_B_y[i]-const.k)/const.k))

# print(const.k)
# # fig, axs = plt.subplots(2, 1, figsize=(6.5,5))
# # axs[0].plot(t[num], savgol_filter(x[num], 51, 3), linewidth=1, label=r'Daten $x$')
# # axs[0].set_xlabel(r'$t$ [s]')
# # axs[0].set_ylabel(r'$x$ [$\mu$m]')
# # axs[0].legend(loc='best')
# # axs[1].plot(t[num], savgol_filter(y[num], 51, 3), linewidth=1, label=r'Daten $y$')
# # axs[1].set_xlabel(r'$t$ [s]')
# # axs[1].set_ylabel(r'$y$ [$\mu$m]')
# # axs[1].legend(loc='best')

# # plt.tight_layout()
# # plt.savefig('plots/pos_keineKraft.pdf')
# # # plt.show()
# # plt.close()

## b)
# current = [50, 70, 100, 150, 200]
current = [40, 50, 70, 100, 200]
print(current)
power = laserleistung(current)
frequenz = []
PSD_x = []
PSD_y = []

for i in np.arange(len(current)):
    with open('data/k_mitKraft_' + str(current[i]) + 'mA_1Hz.FDdat', 'rb') as f:
        clean_lines = (line.replace(b',',b'') for line in f)
        frequenz_temp, PSD_x_temp, PSD_y_temp, egal = np.genfromtxt(clean_lines, unpack=True)
        frequenz.append(frequenz_temp)
        PSD_x.append(PSD_x_temp)
        PSD_y.append(PSD_y_temp)

f_min = [5, 5, 5, 5, 5]
f_max = [50000, 50000, 50000, 50000, 50000]
PSD_x_params = np.zeros((len(current), 2))
PSD_y_params = np.zeros((len(current), 2))
for i in np.arange(len(current)):
    PSD_x_params[i], _ = curve_fit(PSD, frequenz[i][f_min[i]:f_max[i]], PSD_x[i][f_min[i]:f_max[i]])
    PSD_y_params[i], _ = curve_fit(PSD, frequenz[i][f_min[i]:f_max[i]], PSD_y[i][f_min[i]:f_max[i]])


# Fallensteifigkeit:
k_x = computeK(PSD_x_params[:,1])
k_y = computeK(PSD_y_params[:,1])

print('================== Mit Kraft:')
for i in np.arange(len(current)):
    print('x- und y-Roll-Off-Frequenz für ' + str(current[i]) + 'mA: {0:.3e}     {1:.3e}'.format(PSD_x_params[i,1], PSD_y_params[i,1]))

for i in np.arange(len(current)):
    print('x- und y-Fallensteifigkeit für ' + str(current[i]) + 'mA: {0:.3e}     {1:.3e}'.format(k_x[i], k_y[i]))

k_x_params, _ = curve_fit(linear, power, k_x)
k_y_params, _ = curve_fit(linear, power, k_y)

fig, ax = plt.subplots(figsize=(5.5,3.5))
ax.plot(power, k_x, 'x', label=r'Daten $k_x$')
ax.plot(power, k_y, 'x', label=r'Daten $k_y$')
power_fit = np.linspace(0, 100)
ax.plot(power_fit, linear(power_fit, *k_x_params), c='C0', ls='dashed', label='Linearer Fit')
ax.plot(power_fit, linear(power_fit, *k_y_params), c='C1', ls='dashed', label='Linearer Fit')

ax.set_xlabel(r'$P$ [mW]')
ax.set_ylabel(r'$k$ [N/m]')
ax.legend(loc='best')

plt.tight_layout()
plt.savefig('plots/k_mitKraft.pdf')
plt.show()
plt.close()

# Boltzmann-Konstante
t = []
x = []
y = []
for i in np.arange(len(current)):
    t_temp, x_temp, y_temp, egal = np.genfromtxt('data/k_mitKraft_' + str(current[i]) + 'mA_1Hz.TDdat', unpack=True)
    t.append(t_temp)
    x.append(x_temp)
    y.append(y_temp)

x = x / m_x_mean
y = y / m_y_mean

k_B_x = computeK_B(k_x, x*10**(-6))
k_B_y = computeK_B(k_y, y*10**(-6))

for i in np.arange(len(current)):
    print('x- und y-Boltzmannkonst. für ' + str(current[i]) + 'mA: {0:.3e}     {1:.3e}'.format(k_B_x[i], k_B_y[i]))

for i in np.arange(len(current)):
    print('x- und y-Abweichungen für ' + str(current[i]) + 'mA: {0:.3f}     {1:.3f}'.format((k_B_x[i]-const.k)/const.k, (k_B_y[i]-const.k)/const.k))

print(const.k)



