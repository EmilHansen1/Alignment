# %%
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

## ----- Units ----- ##

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

## ----- Functions for constructing the Hamiltonians and Liouvillians ----- ##

def rot_energy(j, B):
    ''' The rotational energy in a.u. '''
    return B * j*(j + 1)


def electric_field(t, e_field_amplitude, fwhm):
    ''' The temporal intensity profile of the Gaussian beam '''
    return e_field_amplitude * np.exp(-two_log_two * (t/fwhm)**2)


def electric_field_squared(t, e_field_amplitude, fwhm):
    ''' The squared electric field '''
    return e_field_amplitude**2 * np.exp(-four_log_two * (t/fwhm)**2)


def alpha(l, m):
    ''' Support function for calcualting cos² matrix elements '''
    return np.sqrt((l - m + 1)*(l + m + 1)/((2*l + 1)*(2*l + 3)))


def beta(l, m):
    ''' Support function for calcualting cos² matrix elements '''
    return np.sqrt((l - m)*(l + m)/((2*l - 1)*(2*l + 1)))


def cos2(lp, mp, l, m):
    if mp == m and lp == l + 2:
        return alpha(l, m)*alpha(l+1, m)
    elif mp == m and lp == l:
        return (alpha(l, m)*beta(l+1, m) + alpha(l-1, m)*beta(l, m))
    elif mp == m and lp == l - 2:
        return beta(l, m)*beta(l-1, m)
    else:
        return 0 


def create_index_dicts(j_max):
    ''' Creates dictionaries for easy indexing '''
    idx_to_lm = {}
    lm_to_idx = {}
    i = 0
    for j in range(j_max + 1):
        for m in range(-j, j + 1):
            idx_to_lm[i] = (j, m)
            lm_to_idx[(j, m)] = i
            i += 1
    
    return idx_to_lm, lm_to_idx


def liouville_norm(rho_vec, dim):
    ''' Calcualtes the norm in Liouville space, Tr[q] '''
    return np.real(np.trace(rho_vec.reshape(dim, dim)))


# %%
## ----- Simulation parameters and settings ----- ##

# Hilbert space settings
j_max = 10
space_dimension = (j_max + 1)**2
idx_to_jm, jm_to_idx = create_index_dicts(j_max)
krylov_n = 15

# Physical parameters (currently for I2 - taken from the alignent calculator)
B_SI = 1.1178 * 2*np.pi * 1e9  # rad Hz 
B = B_SI * au_frequency_per_SI
delta_alpha_SI = 6.0993  # Å³
delta_alpha = delta_alpha_SI * au_per_angstrom**3
I0_SI = 1e12 * 1e4  # W/m²
E0 = np.sqrt(2 * I0_SI / (sc.speed_of_light * sc.epsilon_0)) * au_electric_field_strength_per_SI
FWHM = 300 * au_per_fs
temperature = 1.25  # K
even_abundance = 5.0
odd_abundance = 7.0


## ----- Construct Liouvillan superoperator ----- ##

# First the field-free/time-independent part
L_0 = np.zeros((space_dimension**2, space_dimension**2), dtype=complex)
H_0 = np.zeros((space_dimension, space_dimension)) 
I_mat = np.eye(space_dimension)
for i in range(space_dimension):
    H_0[i, i] = rot_energy(idx_to_jm[i][0], B)
L_0 = -1j * (np.kron(H_0, I_mat) -  np.kron(I_mat, H_0))

# Then the laser-dimer interaction Loiouvillan
L_I = np.zeros((space_dimension**2, space_dimension**2), dtype=complex)
H_I = np.zeros((space_dimension, space_dimension))
for i in range(space_dimension):
    for j in range(space_dimension):
        j1, m1 = idx_to_jm[i]
        j2, m2 = idx_to_jm[j]
        H_I[i, j] = cos2(j1, m1, j2, m2)
L_I = -1j * (np.kron(H_I, I_mat) -  np.kron(I_mat, H_I))

L_tot = L_0 + L_I

# %% 
## ----- Solve the coupled ODE ----- ##

def rho_dot(t, rho):
    ''' The differential equation to be solved '''
    return (L_0 - electric_field_squared(t-300*au_per_fs, E0, FWHM)*delta_alpha*L_I/4) @ rho

# The inital state
single_jm = False

if single_jm:
    # Ground state |00>
    rho_0 = np.zeros((space_dimension, space_dimension), dtype=complex).flatten()
    J, M = 1, 0
    idx = jm_to_idx[(1, 0)]
    rho_0[0] = 1
else:
    # Thermal state
    beta = 1/au_energy_per_SI/(sc.Boltzmann*temperature)
    weights = []
    for i in range(space_dimension):
        j = idx_to_jm[i][1]
        weights.append(even_abundance if j % 2 == 0 else odd_abundance)
    weights = np.diag(np.array(weights))
    rho_0 = np.zeros((space_dimension, space_dimension), dtype=complex)
    rho_0 += sl.expm(-beta*(weights @ H_0))
    rho_0 /= np.trace(rho_0)
    plt.imshow(np.abs(rho_0))
    plt.show()
    rho_0 = rho_0.flatten()

# Solve the ODE

#res = solve_ivp(rho_dot, [0, 500*au_per_ps], rho_0, max_step=1*au_per_ps, method='RK45', rtol=1e-5, atol=1e-5)


# Solve in the Krylov subspace
dt = au_per_ps/20
t_lst = np.linspace(0, 500*au_per_ps, 1000)
res_lst = [rho_0]

for i, t in enumerate(t_lst):
    # First build the Krylov subspace vectors and matrix representation
    q_lst = np.zeros((krylov_n, space_dimension**2), dtype=complex)
    q_lst[0] = rho_0
    L_k = np.zeros((krylov_n, krylov_n), dtype=complex)
    L_t = L_0 + electric_field_squared(t + dt/2 - 300*au_per_fs, E0, FWHM)*delta_alpha*L_I/4
    alpha_lst = np.zeros(krylov_n, dtype=complex)
    beta_lst = np.zeros(krylov_n, dtype=complex) # First is zero!

    # Do the Lanczos iterations!
    for j in range(krylov_n-1):
        alpha_lst[j] = np.dot(np.conjugate(q_lst[j]), L_t @ q_lst[j])
        q_lst[j+1] = L_t @ q_lst[j] #- alpha_lst[j]*q_lst[j] - beta_lst[j]*q_lst[j-1]
        #print(q_lst)
        beta_lst[j] = liouville_norm(q_lst[j+1], space_dimension)
        #print(beta_lst)
        #q_lst[j+i] /= 2
    L_k += np.diag(alpha_lst) + np.diag(beta_lst[1:], -1) + np.diag(beta_lst[1:], 1) 
    
#print(L_k)

## ----- Get expectation values ----- ##
'''
cos2_expval = []
for rho in res.y.T:
    rho = np.reshape(rho, (space_dimension, space_dimension))
    cos2_expval.append(np.real(np.trace(H_I @ rho)))

plt.plot(res.t / au_per_ps, cos2_expval)
plt.show()

plt.imshow(np.abs(res.y.T[-1].reshape(space_dimension, space_dimension))**2)
plt.show()
'''