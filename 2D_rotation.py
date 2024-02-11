# %% ----- IMPORTS AND MATPLOTLIB ASTHETHICS ----- %% #

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.linalg as sl
from scipy.integrate import solve_ivp

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

import matplotlib as mpl
upper = mpl.cm.viridis(np.arange(256))
lower = np.ones((int(256/10),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack((lower, upper))
cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

# %% ----- UNITS ----- %% #

au_velocity_per_SI = 1 / sc.physical_constants['atomic unit of velocity'][0]
speed_of_light_au = 1 / sc.fine_structure
vacuum_permtivity_au = 1/(4*np.pi)
au_electric_field_strength_per_SI = 1/5.14220674763e11
au_length_per_SI = 1/sc.physical_constants['atomic unit of length'][0]
au_time_per_SI = 1/sc.physical_constants['atomic unit of time'][0]
au_frequency_per_SI = sc.physical_constants['atomic unit of time'][0]
au_per_angstrom = 1.8897259886
au_per_fs = 1e-15 * au_time_per_SI
au_per_ps = 1e-12 * au_time_per_SI
two_log_two = 2 * np.log(2)
four_log_two = 4 * np.log(2)
au_energy_per_SI = 1 / 4.3597447222071e-18

# %% ----- UTILITY FUNCTIONS AND CLASSES ----- %% #

def rot_energy(B, j):
    ''' The rotational energy for a 2D rotor '''
    return B * j**2


def electric_field(t, e_field_amplitude, fwhm):
    ''' The temporal intensity profile of the Gaussian beam '''
    return e_field_amplitude * np.exp(-two_log_two * (t/fwhm)**2)


def electric_field_squared(t, e_field_amplitude, fwhm):
    ''' The squared electric field '''
    return e_field_amplitude**2 * np.exp(-four_log_two * (t/fwhm)**2)


def cos2_element(j1, j2):
    if np.abs(j1 - j2) == 2:
        return 1/4
    elif np.abs(j1 - j2) == 0:
        return 1/2
    else:
        return 0

# %% ----- PULSE AND TARGET SETTINGS SETTINGS ----- %% #

B_SI = 1.33 * 2*np.pi * 1e9  # rad Hz 1.1178 
B = B_SI * au_frequency_per_SI
delta_alpha_SI = 31.2#6.0993  # Å³
delta_alpha = delta_alpha_SI * au_per_angstrom**3
I0_SI = 2e12 * 1e4  # W/m²
E0 = np.sqrt(2 * I0_SI / (sc.speed_of_light * sc.epsilon_0)) * au_electric_field_strength_per_SI
FWHM = 300 * au_per_fs
T = 0

# %% ----- CREATE MATRICES ----- %% #

# Number of |j> basis states
j_max = 50
n_j = 2*j_max + 1

# Create interaction Hamiltonian
j_lst = np.arange(-j_max, j_max + 1)
cos2 = np.zeros((n_j, n_j))
H_I = np.zeros((n_j, n_j))
H_0 = np.zeros((n_j, n_j))

for i, j1 in enumerate(j_lst):
    H_0[i,i] = rot_energy(B, j1)
    for j, j2 in enumerate(j_lst):
        cos2[i,j] = cos2_element(j1, j2)
        H_I[i,j] = -0.25 * delta_alpha * cos2[i,j]


def get_cos2(res):
    cos2_expval = []
    for rho_arr in res.y.T:
        rho = np.reshape(rho_arr, (n_j, n_j))
        cos2_expval.append(np.real(np.trace(cos2 @ rho)))
    return cos2_expval

# Create ODE for LvN equation
def solve_LvN_eqn(T, t_end, fwhm, ratio=0, delta=0):

    # Create inital state
    rho_0 = np.zeros((n_j, n_j), dtype=complex)
    if T == 0:
        rho_0[j_max, j_max] = 1
    else:
        beta = 1/au_energy_per_SI/(sc.Boltzmann*T)
        rho_0 += sl.expm(-beta*H_0)
        rho_0 /= np.trace(rho_0)
        #plt.imshow(np.abs(rho_0))

    def rho_dot(t, rho):
        rho = rho.reshape((n_j, n_j))
        H = H_0 + H_I * (electric_field_squared(t, E0, fwhm))
        return (-1j * (H @ rho - rho @ H)).flatten()
    
    res = solve_ivp(rho_dot, [-2*fwhm, t_end], rho_0.flatten(), max_step=5*au_per_ps, method='RK45', rtol=1e-7, atol=1e-7)
    return res



# %% ----- SOLVE THE LvN EQUATION ----- %% #

T = 0
res = solve_LvN_eqn(T, 200*au_per_ps, 300*au_per_fs)
cos2_trace = get_cos2(res)
plt.plot(res.t, cos2_trace)

# %% ----- PLOT ALIGNMENT ----- %% #

echo_y = []
echo_x = []
diff = []
for i, (res, delta) in enumerate(zip(results, deltas)):
    #plt.plot(res.t / au_per_ps, get_cos2(res))
    mask = (res.t < (2*delta + 5)*au_per_ps) & ((2*delta - 5)*au_per_ps < res.t)
    echo_x.append((res.t[mask] - 2*delta*au_per_ps) / au_per_ps)
    echo_y.append(np.array(get_cos2(res))[mask])
    diff.append(np.max(echo_y[i]) - np.min(echo_y[i]))
    plt.plot(echo_x[i], echo_y[i])
plt.show()

# %%
for i in range(14):
    plt.plot(echo_x[i], echo_y[i])

# %%

plt.plot(deltas[:8], diff[:8], 'ko')


# %% ----- PLOT THE DENSITY MATRIX ----- %% #

plt.imshow(np.abs(np.reshape(res.y.T[-1], (n_j, n_j)))**2, cmap=cmap)
#plt.xlim(20, 80)
#plt.ylim(20, 80)
np.trace(np.reshape(res.y.T[-1], (n_j, n_j)))

# %% ----- PLOT POPULATIONS ----- %% #

plt.plot(res.t / au_per_ps, [np.reshape(rho, (n_j, n_j))[j_max, j_max] for rho in res.y.T])
plt.plot(res.t / au_per_ps, [np.reshape(rho, (n_j, n_j))[j_max+2, j_max+2] for rho in res.y.T])
plt.plot(res.t / au_per_ps, [np.reshape(rho, (n_j, n_j))[j_max+4, j_max+4] for rho in res.y.T])
plt.plot(res.t / au_per_ps, [np.reshape(rho, (n_j, n_j))[j_max+6, j_max+6] for rho in res.y.T])
plt.plot(res.t / au_per_ps, [np.reshape(rho, (n_j, n_j))[j_max+8, j_max+8] for rho in res.y.T])
plt.plot(res.t / au_per_ps, [np.reshape(rho, (n_j, n_j))[j_max+10, j_max+10] for rho in res.y.T])
plt.xlim(-1, 1)

# %%
