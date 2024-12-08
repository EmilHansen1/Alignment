import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from scipy.fft import fft

# Matplotlib asthetics
major = 8
minor = 4
width = 0.75
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


fig, ax = plt.subplots(1, 1)
times, cos2, field = np.loadtxt('cos2.out', delimiter=',').T
#cos2 = (cos2 - 0.5) * 0.68 + 0.5 # For Rb
#cos2 = (cos2 - 0.5) * 0.42 + 0.5 # For K
ax.plot(times, cos2, 'b-')

times, cos2, field = np.loadtxt('cos2_K_singlet.out', delimiter=',').T
ax.plot(times, cos2, 'k--')

ax2 = ax.twinx()
ax2.fill_between(times, field, alpha=0.4, color='gray')
ax2.plot(times, field / np.max(field), c='gray')
ax2.set(ylim=(0, 1.2))
#times, cos2, field = np.loadtxt('cos2_circ.out', delimiter=',').T
#ax.plot(times, cos2, 'k--')

#times, cos2 = np.loadtxt('Data/Rb2_singlet_trace.csv', delimiter=',').T
#ax.plot(times, cos2 + (0.5 - 0.465369), 'k-')


plt.show()

#times, cos2_1, field = np.loadtxt('cos2_Rb_1.out', delimiter=',').T
#times, cos2_2, field = np.loadtxt('cos2_Rb_2.out', delimiter=',').T
#plt.plot(times, (cos2_1*72.2 + cos2_2*27.8) / (93.2581 + 6.7302), 'b-', label='My 2D code')

#plt.plot(times, cos2_1)
#plt.plot(times, cos2_2)


#iters, times, cos2 = np.loadtxt(r'Tests/K3DWithUInhomogenousB0.66D0.02I2.59e10WFA.csv', delimiter=',', skiprows=1).T
#iters, times, cos2 = np.loadtxt(r'Tests/Rb3DWithUInhomogenousB0.28D0.02I2.81e10WFA.csv', delimiter=',', skiprows=1).T
#iters, times, cos2 = np.loadtxt(r'Tests/K2DB0.69I2.59e10.csv', delimiter=',', skiprows=1).T
#iters, times, cos2 = np.loadtxt(r'Tests/Rb2DB0.28I2.81e10.csv', delimiter=',', skiprows=1).T
#iters, times, cos2 = np.loadtxt(r'Tests/Na3DSingletWithUInhomogenousB4.53D0.19I3.53e10WFA.csv', delimiter=',', skiprows=1).T
#iters, times, cos2 = np.loadtxt(r'Tests/NaSingletTest.csv', delimiter=',', skiprows=1).T
'''
plt.figure(figsize=(8, 3))
plt.plot(times, cos2, 'k-', label='Areg\'s 2D code')
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle \text{cos}^2\theta_\text{2D}\rangle$')
plt.minorticks_on()
#plt.legend(frameon=False)
#plt.savefig('Rb 2D comparison.png', dpi=300)
plt.show()
'''
'''ax2 = plt.gca().twinx()
ax2.plot(times, field, 'gray')
ax2.fill_between(times, field, color='gray', alpha=0.5)
ax2.set_ylabel(r'$I(t)/I_0$')
ax2.set_ylim(0, 1.3)
ax2.minorticks_on()'''


'''
fig, axs = plt.subplots(3, 1, sharex='all', figsize=(5, 8), )

# Na2
times_Na, cos2_Na, _ = np.loadtxt('cos2_Na.out', delimiter=',').T
axs[0].plot(times_Na, cos2_Na, 'b-', lw=0.75)

# K2
iters, times, cos2 = np.loadtxt(r'Tests/K2DB0.69I2.59e10.csv', delimiter=',', skiprows=1).T
axs[1].plot(times, cos2, 'k-', lw=1.5)
times_K, cos2_K_1, _ = np.loadtxt('cos2_K_1.out', delimiter=',').T
times_K, cos2_K_2, _ = np.loadtxt('cos2_K_2.out', delimiter=',').T
p1 = 0.932581
p2 = 0.067302
axs[1].plot(times_K, cos2_K_1* p1 + cos2_K_2 * p2, 'b-', lw=0.75)

# Rb2
iters, times, cos2 = np.loadtxt(r'Tests/Rb2DB0.28I2.81e10.csv', delimiter=',', skiprows=1).T
axs[2].plot(times, cos2, 'k-', lw=1.5)
times_Rb, cos2_Rb_1, _ = np.loadtxt('cos2_Rb_1.out', delimiter=',').T
times_Rb, cos2_Rb_2, _ = np.loadtxt('cos2_Rb_2.out', delimiter=',').T
p1 = 0.722
p2 = 0.278
axs[2].plot(times_Rb, cos2_Rb_1 * p1 + cos2_Rb_2 * p2, 'b-', lw=0.75)

labels = [r'Na$_2\,(^3\Sigma_\text{u}^+)$', r'K$_2\,(^3\Sigma_\text{u}^+)$', r'Rb$_2\,(^3\Sigma_\text{u}^+)$']
for ax, label in zip(axs, labels):
    ax.minorticks_on()
    ax.legend(frameon=False)
    ax.set(xlabel = None if ax is not axs[2] else r'$t$ (ps)', ylabel=r'$\langle\cos^2\theta_\text{2D}\rangle$')
    t = ax.text(0.95, 0.9, label, transform=ax.transAxes, fontsize=16, horizontalalignment='right', verticalalignment='top')
    t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='black'))


plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('Ak2_triplets.png', dpi=300)
plt.show()
'''

'''
if __name__ == '__main__':
    fig, axs = plt.subplots(3, 1, figsize=(5, 8), sharex='all', constrained_layout=True)
    for i, (ax, fname) in enumerate(zip(axs, ['Li', 'Na', 'K'])):
        times, cos2, field = np.loadtxt(f'cos2_{fname}.out', delimiter=',').T
        ax.plot(times, cos2)
        ax.minorticks_on()
        ax.set_ylabel(r'$\langle \text{cos}^2\theta_\text{2D}\rangle$')

        ax2 = ax.twinx()
        ax2.fill_between(times, field, color='gray', alpha=0.5)
        ax2.fill_between(times, field, color='gray', alpha=0.5)
        ax2.set_ylabel(r'$I(t)/I_0$')
        ax2.set_ylim(0, 1.3)
    axs[-1].set_xlabel('Time (ps)')
    plt.savefig('triplet_simulations.png', dpi=500)
    plt.show()
'''