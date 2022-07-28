from random import gauss
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as const

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern Roman"],
    "font.size": 7.0})

def gaussFit(t, A, delta_t):
    return A * np.exp(-4*np.log(2)*(t/delta_t)**2)

def gaussFit1(t, A, delta_t, t_0):
    return A * np.exp(-4*np.log(2)*((t-t_0)/delta_t)**2)


################################################################################################
## Original
lam, I_spek = np.genfromtxt('data/spektrum_original.SSM', skip_header=3, unpack=True)
t, I_auto = np.genfromtxt('data/autokorellation_original.txt', usecols=(1,2), unpack=True)
I_auto = I_auto/I_auto.max()
i_max = np.argmax(I_auto)
t = t-t[i_max]   # in mm
t = 2 * t / (1000 * const.c) * 10**(12)   # in ps

fig, axs = plt.subplots(1, 2, figsize=(6, 2.8))
axs[0].plot(lam, I_spek/I_spek.max(), label='Spektrum')
axs[0].axvline(1550, ls='dashed', c='grey', lw=1.5, alpha=0.9, label='Zentrale Wellenlänge\n'+r'$1550\,$nm')
axs[0].set_xlabel(r'$\lambda$ [nm]')
axs[0].set_ylabel(r'$I$ [a.u.]')
axs[0].legend(bbox_to_anchor=(-0.01, 1.3), loc='upper left')

i_fit = 15
params, cov = curve_fit(gaussFit, t[i_max-i_fit:i_max+i_fit], I_auto[i_max-i_fit:i_max+i_fit])
print('Halbwertsbreite original:  {0:.1f}fs'.format(np.abs(params[1])*1000 /np.sqrt(2)))
print(np.abs(np.abs(params[1])*1000 /np.sqrt(2) - 100)/100)

axs[1].plot(t, I_auto, label='Autokorrelationsspur')
axs[1].plot(t, gaussFit(t, *params), 'r--', lw=1, label='Gauss-Fit')
axs[1].plot(t[i_max-i_fit:i_max+i_fit], I_auto[i_max-i_fit:i_max+i_fit], c='#1f77b4', lw=6, alpha=0.4, label='gefittete Werte')
axs[1].set_xlabel(r'$t$ [ps]')
axs[1].set_ylabel(r'$I$ [a.u.]')
axs[1].set_xlim(-0.7, 0.7)
axs[1].legend(bbox_to_anchor=(-0.01, 1.32), loc='upper left')
plt.tight_layout()
# plt.subplots_adjust(wspace=0.2)
plt.savefig('plots/Original.pdf')
# plt.show()
plt.close()


################################################################################################
## Mit 30nm-Bandpassfilter
lam, I_spek = np.genfromtxt('data/spektrum_bandpass_1550_30.SSM', skip_header=3, unpack=True)
t, I_auto = np.genfromtxt('data/autokorrelation_bandpass_1550_30.txt', usecols=(1,2), unpack=True)
I_auto = I_auto/I_auto.max()
i_max = np.argmax(I_auto)
t = t-t[i_max]   # in mm
t = 2 * t / (1000 * const.c) * 10**(12)   # in ps

fig, axs = plt.subplots(1, 2, figsize=(6, 2.8))
axs[0].plot(lam, I_spek/I_spek.max(), label='Spektrum')
axs[0].axvline(1550, ls='dashed', c='grey', lw=1.5, alpha=0.9, label='Zentrale Wellenlänge\n'+r'$1550\,$nm')
axs[0].set_xlabel(r'$\lambda$ [nm]')
axs[0].set_ylabel(r'$I$ [a.u.]')
axs[0].legend(bbox_to_anchor=(-0.01, 1.3), loc='upper left')

i_fit = 15
params, cov = curve_fit(gaussFit, t[i_max-i_fit:i_max+i_fit], I_auto[i_max-i_fit:i_max+i_fit])
print('Halbwertsbreite 30nm-Bandpass:  {0:.1f}fs'.format(np.abs(params[1])*1000 /np.sqrt(2)))

axs[1].plot(t, I_auto, label='Autokorrelationsspur')
axs[1].plot(t, gaussFit(t, *params), 'r--', lw=1, label='Gauss-Fit')
axs[1].plot(t[i_max-i_fit:i_max+i_fit], I_auto[i_max-i_fit:i_max+i_fit], c='#1f77b4', lw=6, alpha=0.4, label='gefittete Werte')
axs[1].set_xlabel(r'$t$ [ps]')
axs[1].set_ylabel(r'$I$ [a.u.]')
axs[1].set_xlim(-1, 1)
axs[1].legend(bbox_to_anchor=(-0.01, 1.32), loc='upper left')
plt.tight_layout()
plt.subplots_adjust(wspace=0.2)
plt.savefig('plots/30nm_Bandpass.pdf')
# plt.show()
plt.close()


################################################################################################
## Mit 12nm-Bandpassfilter
lam, I_spek = np.genfromtxt('data/spektrum_bandpass_1550_12.SSM', skip_header=3, unpack=True)
t, I_auto = np.genfromtxt('data/autokorrelation_bandpass_1550_12.txt', usecols=(1,2), unpack=True)
I_auto = I_auto/I_auto.max()
i_max = np.argmax(I_auto)
t = t-t[i_max]   # in mm
t = 2 * t / (1000 * const.c) * 10**(12)   # in ps

fig, axs = plt.subplots(1, 2, figsize=(6, 2.8))
axs[0].plot(lam, I_spek/I_spek.max(), label='Spektrum')
axs[0].axvline(1550, ls='dashed', c='grey', lw=1.5, alpha=0.9, label='Zentrale Wellenlänge\n'+r'$1550\,$nm')
axs[0].set_xlabel(r'$\lambda$ [nm]')
axs[0].set_ylabel(r'$I$ [a.u.]')
axs[0].legend(bbox_to_anchor=(-0.01, 1.3), loc='upper left')

i_fit = 15
params, cov = curve_fit(gaussFit, t[i_max-i_fit:i_max+i_fit], I_auto[i_max-i_fit:i_max+i_fit])
print('Halbwertsbreite 12nm-Bandpass:  {0:.1f}fs'.format(np.abs(params[1])*1000 /np.sqrt(2)))

axs[1].plot(t, I_auto, label='Autokorrelationsspur')
axs[1].plot(t, gaussFit(t, *params), 'r--', lw=1, label='Gauss-Fit')
axs[1].plot(t[i_max-i_fit:i_max+i_fit], I_auto[i_max-i_fit:i_max+i_fit], c='#1f77b4', lw=6, alpha=0.4, label='gefittete Werte')
axs[1].set_xlabel(r'$t$ [ps]')
axs[1].set_ylabel(r'$I$ [a.u.]')
axs[1].set_xlim(-1, 1)
axs[1].legend(bbox_to_anchor=(-0.01, 1.32), loc='upper left')
plt.tight_layout()
plt.subplots_adjust(wspace=0.2)
plt.savefig('plots/12nm_Bandpass.pdf')
# plt.show()
plt.close()


################################################################################################
## Mit Si-Block
lam, I_spek = np.genfromtxt('data/spektrum_Si_12mm.SSM', skip_header=3, unpack=True)
t, I_auto = np.genfromtxt('data/autokorrelation_Si_12mm.txt', usecols=(1,2), unpack=True)
I_auto = I_auto/I_auto.max()
i_max = np.argmax(I_auto)
t = t-t[i_max]   # in mm
t = 2 * t / (1000 * const.c) * 10**(12)   # in ps

fig, axs = plt.subplots(1, 2, figsize=(6, 2.8))
axs[0].plot(lam, I_spek/I_spek.max(), label='Spektrum')
axs[0].axvline(1550, ls='dashed', c='grey', lw=1.5, alpha=0.9, label='Zentrale Wellenlänge\n'+r'$1550\,$nm')
axs[0].set_xlabel(r'$\lambda$ [nm]')
axs[0].set_ylabel(r'$I$ [a.u.]')
axs[0].legend(bbox_to_anchor=(-0.01, 1.3), loc='upper left')

i_fit = 25
i_schief = 0
params, cov = curve_fit(gaussFit, t[i_max-i_fit+i_schief:i_max+i_fit+i_schief], I_auto[i_max-i_fit+i_schief:i_max+i_fit+i_schief])
# params_u = unp.uarray(params, np.sqrt(np.diag(cov)))
print('Halbwertsbreite 12mm-Si-Block:  {0:.1f}fs'.format(np.abs(params[1])*1000 /np.sqrt(2)))

axs[1].plot(t, I_auto, label='Autokorrelationsspur')
axs[1].plot(t, gaussFit(t, *params), 'r--', lw=1, label='Gauss-Fit')
axs[1].plot(t[i_max-i_fit+i_schief:i_max+i_fit+i_schief], I_auto[i_max-i_fit+i_schief:i_max+i_fit+i_schief], c='#1f77b4', lw=6, alpha=0.4, label='gefittete Werte')
axs[1].set_xlabel(r'$t$ [ps]')
axs[1].set_ylabel(r'$I$ [a.u.]')
axs[1].set_xlim(-2.2, 2.2)
axs[1].legend(bbox_to_anchor=(-0.01, 1.32), loc='upper left')
plt.tight_layout()
plt.subplots_adjust(wspace=0.2)
plt.savefig('plots/Si_12mm.pdf')
# plt.show()
plt.close()


################################################################################################
## Mit SF11 Glas-Blöcken
name_numbers1 = [2, 5, 9, 13, 17, 23, 29]
name_numbers2 = [85, 31, 63, 51, 64, 85, 83]

t = []
I_auto = []
for i in np.arange(len(name_numbers1)):
    t_temp, I_auto_temp = np.genfromtxt('data/autokorrelation_Glas_'+str(name_numbers1[i])+'_'+str(name_numbers2[i])+'mm.txt', usecols=(1,2), unpack=True)
    t.append(t_temp)
    I_auto.append(I_auto_temp)

i_max = np.empty(len(name_numbers1))
for i in np.arange(len(name_numbers1)):
    I_auto[i] = I_auto[i]/I_auto[i].max()
    i_max[i] = np.argmax(I_auto[i])
    t[i] = t[i]-t[i][int(i_max[i])]   # in mm
    t[i] = 2 * t[i] / (1000 * const.c) * 10**(12)   # in ps

t_fit = np.linspace(-0.2, 0.2, 1000)
fig, ax = plt.subplots(figsize=(6, 3.5))
for i in np.arange(len(name_numbers1)):
    i_fit = 100
    params, cov = curve_fit(gaussFit1, t[i][int(i_max[i])-i_fit:int(i_max[i])+i_fit], I_auto[i][int(i_max[i])-i_fit:int(i_max[i])+i_fit])
    params_u = unp.uarray(params, np.sqrt(np.diag(cov)))
    # ax.plot(t[i], I_auto[i], lw=0.5, label='{0},{1}$\,$mm'.format(name_numbers1[i], name_numbers2[i]))
    ax.plot(t_fit, gaussFit1(t_fit+params[2], *params), lw=0.5, label='{0},{1}$\,$mm'.format(name_numbers1[i], name_numbers2[i]))
    print('Halbwertsbreite {0},{1}mm-Glass:  {2:.1f}fs'.format(name_numbers1[i], name_numbers2[i], np.abs(params_u[1])*1000 /np.sqrt(2)))

ax.set_xlabel(r'$t$ [ps]')
ax.set_ylabel(r'$I$ [a.u.]')
ax.set_xlim(-0.2, 0.2)
ax.legend(loc='upper left')
plt.tight_layout()
plt.savefig('plots/Glas.pdf')
# plt.show()
plt.close()




