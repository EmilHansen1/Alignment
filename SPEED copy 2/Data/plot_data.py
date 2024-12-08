import matplotlib.pyplot as plt
import numpy as np
import matplotlib

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

# It's plots o'clock
fig, axs = plt.subplots(2, 1, constrained_layout=True)

# Plot K2 stuffs
# B_2D = 1.78 GHz
times, cos2_K2_1, _ = np.loadtxt('cos2_K2_1.out', delimiter=',').T
times, cos2_K2_2, _ = np.loadtxt('cos2_K2_2.out', delimiter=',').T
exp_times, exp_cos2 = np.loadtxt('K2_singlet_trace.csv', delimiter=',').T
abundance_1, abundance_2 = 0.933, 0.087
theo_cos2 = (cos2_K2_1 * abundance_1 + cos2_K2_2 * abundance_2 - 0.5) * 0.38 + 0.5
axs[0].plot(times, theo_cos2, 'b-')
axs[0].plot(exp_times, exp_cos2, 'k-')

# Plot Rb2 stuffs
# B_2D = 0.69 GHz
times, cos2_Rb2_1, _ = np.loadtxt('cos2_Rb2_1.out', delimiter=',').T
times, cos2_Rb2_2, _ = np.loadtxt('cos2_Rb2_2.out', delimiter=',').T
exp_times, exp_cos2 = np.loadtxt('Rb2_singlet_trace.csv', delimiter=',').T
abundance_1, abundance_2 = 0.722, 0.278
theo_cos2 = (cos2_Rb2_1 * abundance_1 + cos2_Rb2_2 * abundance_2 - 0.5) * 0.7 + 0.5
axs[1].plot(times, theo_cos2, 'b-')
axs[1].plot(exp_times, exp_cos2 + (0.5 - 0.465369), 'k-')

for ax in axs:
    ax.minorticks_on()
    ax.set(xlabel=r'$t$ (ps)' if ax is axs[1] else None, ylabel=r'$\langle\cos^2\theta_\text{2D}\rangle$')
plt.savefig('plot.png', dpi=350)
plt.show()



