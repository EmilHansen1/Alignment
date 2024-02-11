import matplotlib.pyplot as plt
import numpy as np

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

time, cos2 = np.loadtxt('cos2.out').T
plt.plot(time, cos2)
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle \text{cos}^2\theta_\text{2D}\rangle$')
plt.minorticks_on()
plt.show()