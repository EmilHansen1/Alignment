# %% Imports

import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.constants as sc

# Matplotlib asthetics
major = 6
minor = 3
width = 1
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
plt.rc("axes", labelsize=16) # 18
plt.rc("xtick", labelsize=14, top=True, direction="in")
plt.rc("ytick", labelsize=14, right=True, direction="in")
plt.rc("axes", titlesize=18)
plt.rc("legend", fontsize=14)
plt.rcParams['font.family'] = "serif"
plt.rcParams['axes.linewidth'] = width
plt.rcParams['xtick.minor.width'] = width
plt.rcParams['xtick.major.width'] = width
plt.rcParams['ytick.minor.width'] = width
plt.rcParams['ytick.major.width'] = width
plt.rcParams['xtick.major.size'] = major
plt.rcParams['xtick.minor.size'] = minor
plt.rcParams['ytick.major.size'] = major
plt.rcParams['ytick.minor.size'] = minor

# Conversion factors
au_per_s = 2.4188843265857e17 # 4.13413732e16 
au_per_ps = 1e-12*au_per_s 


# %% Parameters 

ang_per_au = 0.529177249
au_per_s = 4.13413732e16

B = 1.1178*1e9 / au_per_s #1.700125755756e-7 / (2*np.pi) # 0.4600310255208 in GHz
I = 3.0e13
j = 0
m = 0
j_max = 40
a_parr = 14.562 #/ ang_per_au**3 #98.05
a_perp = 8.4627 #/ ang_per_au**3 #52.779
fwhm = 100
offset = 20
ratio = 1


cmd = f'build/rot_tdse_solver {B} {I} {j} {m} {j_max} {a_parr} {a_perp} {fwhm} {offset} {ratio}'

os.system(cmd)

t_list, cos2_list = np.loadtxt('dat.out').T

t_list = t_list/au_per_ps

plt.plot(t_list, cos2_list, color='darkslateblue')
plt.minorticks_on()
plt.xlabel(r'Delay $\tau$ (ps)')
plt.ylabel(r'$\langle\cos^2\theta\rangle$')
plt.xlim(0, 500)
plt.show()


'''
# %% Boltzman ensemble
temperature = 0.10
B_SI = 1.11863*1e9

partition_function = np.sum([(2*J + 1) * np.exp(-B_SI*J*(J + 1)*sc.hbar/(sc.Boltzmann * temperature)) for J in range(0, j_max + 1)])
weights = np.array([np.exp(-B_SI*J*(J + 1)*sc.hbar/(sc.Boltzmann * temperature)) for J in range(0, j_max + 1)]) / partition_function

# %%

cos2_thermal = []
n_sims = (j_max + 1)**2
counter = 0
for i in range(0, 7):
    w = weights[i]
    for k in range(-i, i+1):
        cmd = f'./Alignment {B} {I} {i} {k} {j_max} {a_parr} {a_perp} {fwhm} {offset} {ratio}'
        os.system(cmd)
        t_list, cos2 = np.loadtxt('dat.out').T
        t_list = t_list/au_per_ps        
        cos2_thermal.append(w*cos2)
        plt.plot(t_list, cos2)
        counter += 1
        print(f'Sim. {counter}/{n_sims} done!')
cos2_thermal = np.array(cos2_thermal)

print('Done!')


# %%

plt.plot(t_list, cos2_thermal.T)
plt.show()

# %%

plt.plot(t_list, np.sum(cos2_thermal, axis=0))
plt.xlim(0, 300)
plt.ylim(0.47, 0.5)
plt.xlabel(r'Delay $\tau$ (ps)')
plt.ylabel(r'$\langle\cos^2\theta\rangle$')
plt.minorticks_on()
plt.show()
# %%
plt.plot(weights, 'ko', mfc='w')

# %%
'''