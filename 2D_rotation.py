# %% ----- IMPORTS AND MATPLOTLIB ASTHETHICS ----- %% #

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.linalg as sl
from scipy.integrate import solve_ivp
import scipy.sparse as sp
from scipy.signal import savgol_filter
from scipy.interpolate import CubicSpline
import time

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
vacuum_permtivity_au = 1 / (4*np.pi)
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


def peak_intensity(E, tau, wx, wy):
    return 4 * np.sqrt(np.log(2)) * E / (np.pi**(3/2) * tau * wx * wy)

print(peak_intensity(170e-6, 175e-12, 80.54*1e-6/2,  179.4*1e-6/2)*1e-4)


def cos2_element(j1, j2):
    if np.abs(j1 - j2) == 2:
        return 1/4
    elif np.abs(j1 - j2) == 0:
        return 1/2
    else:
        return 0
    

def multiply_cos2_2D(psi): 
    ''' Quick matrix multiplication with tri-diagonal matrix '''
    dim = len(psi)
    result = np.zeros(dim, dtype=complex)
    result[0] += psi[0]/2.0 + psi[2]/4.0
    result[1] += psi[1]/2.0 + psi[3]/4.0
    for i in range(2, dim - 1):
        result[i] += psi[i-2]/4.0 + psi[i]/2.0 + psi[i+2]/4.0
    
    result[dim-2] = psi[dim-4]/2.0 + psi[dim-2]/4.0
    result[dim-1] = psi[dim-3]/2.0 + psi[dim-1]/4.0
    return result
    


# %% ----- PULSE AND TARGET SETTINGS SETTINGS ----- %% #

B_SI = 1.33 * 2*np.pi * 1e9  # rad Hz 1.1178 
B = B_SI * au_frequency_per_SI
delta_alpha_SI = 31.2 #6.0993  # Å³
delta_alpha = delta_alpha_SI * au_per_angstrom**3
I0_SI = 8e9 * 1e4  # W/m²
E0 = np.sqrt(2 * I0_SI / (sc.speed_of_light * sc.epsilon_0)) * au_electric_field_strength_per_SI
FWHM = 185 * au_per_ps
T = 0

# %% ----- CREATE MATRICES ----- %% #

'''# Number of |j> basis states
j_max = 30
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

T = 0.4
res = solve_LvN_eqn(T, 200*au_per_ps, 300*au_per_fs)
cos2_trace = get_cos2(res)
plt.plot(res.t / au_per_ps, cos2_trace)'''


# %% ----- CREATE MATRICES ----- %% #

# Load pulse
x, pulse, _, _ = np.loadtxt(r'Pulses/Li_26_8_24.csv', delimiter=',', skiprows=1).T
pulse -= min(pulse)
pulse /= max(pulse)
x0 = 25.192402405721207 - 2
pulse_time = (x - x0) / sc.speed_of_light  / 1e-12 / 1000 * 4
pulse_lims = (-220, 250)
#pulse = savgol_filter(pulse, 2, 1)
plt.plot(pulse_time, pulse)

def xc_pulse(t):
    return np.interp(t, pulse_time, pulse) * np.heaviside(t - pulse_lims[0], 0) * np.heaviside(-t + pulse_lims[1], 0)

x_interp = np.linspace(-250, 210, 500)
x_interp_long = np.linspace(-300, 500, 1000)
field = CubicSpline(x_interp_long, savgol_filter(xc_pulse(x_interp_long), 35, 3), )

plt.plot(x_interp, savgol_filter(xc_pulse(x_interp), 35, 3))

x_interp_long = np.linspace(-500, 1000, 1000)
plt.plot(x_interp_long, field(x_interp_long))


# %%

# Number of |j> basis states
j_max = 9#9
n_j = 2*j_max + 1

# Create interaction Hamiltonian
j_lst = np.arange(n_j) - j_max
cos2 = np.zeros((n_j, n_j))
H_I = np.zeros((n_j, n_j))
H_0 = np.zeros((n_j, n_j))

for i, j1 in enumerate(j_lst):
    H_0[i,i] = rot_energy(1, j1)
    for j, j2 in enumerate(j_lst):
        cos2[i,j] = cos2_element(j1, j2)
        H_I[i,j] = -0.25 * delta_alpha * cos2[i,j]


# Sparse versions
cos2_diagonals = [[0.25]*(2*j_max - 1), [0.5]*(2*j_max + 1), [0.25]*(2*j_max - 1)]
cos2_s = sp.diags(cos2_diagonals, [-2, 0, 2], shape=(n_j, n_j), format='dia', dtype=complex)
H_I_s = -0.25 * delta_alpha * cos2_s
H_0_s = sp.diags(np.diag(H_0), 0, shape=(n_j, n_j), format='dia', dtype=complex)


n_BE = lambda w: 1/(np.exp(w/au_energy_per_SI/(sc.Boltzmann*0.38)) - 1)
# Create Lindblad operators and dissipators
# Currently just projection operators, projecting one energy level down
L_ops = []
j = -j_max
for i in range(n_j): # Loop over ALL quantum numbers
    if j != 0:
        L_op = np.zeros((n_j, n_j))
        if j > 0:
            L_op[i-1, i] = 1#np.sqrt(n_BE(np.abs(rot_energy(B, j) - rot_energy(B, j-1))) + 1)
        else:
            L_op[i+1, i] = 1#np.sqrt(n_BE(np.abs(rot_energy(B, j+1) - rot_energy(B, j))) + 1)
        L_ops.append(L_op)
    j += 1    

L_op = np.sum(L_ops, axis=0)
L_op_dag = np.conjugate(L_op.T)

''''
j = -j_max
for i in range(n_j): # Loop over ALL quantum numbers
    if j != 0:
        L_op = np.zeros((n_j, n_j))
        if j > 0:
            L_op[i, i-1] =0.1# np.sqrt(n_BE(np.abs(rot_energy(B, j) - rot_energy(B, j-1))))
        else:
            L_op[i, i+1] = 0.1 #np.sqrt(n_BE(np.abs(rot_energy(B, j+1) - rot_energy(B, j))))
        L_ops.append(L_op)
    j += 1
'''

K_ops = []
for i in range(n_j):
    K_op = np.zeros((n_j, n_j))
    K_op[i, i] = 1
    K_ops.append(K_op) 

def get_cos2(res):
    cos2_expval = []
    for rho_arr in res.y.T:
        rho = np.reshape(rho_arr, (n_j, n_j))
        cos2_expval.append(np.real(np.trace(cos2 @ rho)))
    return np.array(cos2_expval)

def get_cos2_exp(res):
    cos2_expval = []
    for psi in res.y.T:
        cos2_expval.append(np.real(np.dot(np.conj(psi), cos2 @ psi)))
    return cos2_expval


def solve_ode(t_end):
    psi_0 = np.zeros(n_j, dtype=complex)
    psi_0[j_max-1] = 1.0 + 0.0j
    
    def psi_dot(t, psi):
        H = H_0 + H_I * electric_field_squared(t, E0, FWHM)
        return -1j * H @ psi
    
    return solve_ivp(psi_dot, [-2*FWHM, t_end], psi_0, max_step=25*au_per_ps, method='RK45', rtol=1e-7, atol=1e-7)

# %%

B_SI = 5.25 * 2*np.pi * 1e9  # rad Hz 1.1178 1.33
B_SI = 1.29 * 2*np.pi * 1e9
B = B_SI * au_frequency_per_SI
delta_alpha_SI = 31.2#65.73#31.2 #6.0993  # Å³
delta_alpha = delta_alpha_SI * au_per_angstrom**3
I0_SI = 5.86e9 * 1e4  # W/m²
E0 = np.sqrt(2 * I0_SI / (sc.speed_of_light * sc.epsilon_0)) * au_electric_field_strength_per_SI
FWHM = 173 * au_per_ps#1.2 * au_per_ps
n_even = 5
n_odd = 3

gamma = 1/(350 * au_per_ps)
pd = 1/(300 * au_per_ps)

def solve_rho(t_end, temp, I0, B, plot=False):
    # Calculate intensity
    I0_SI = I0 * 1e4  # W/m²
    E0 = np.sqrt(2 * I0_SI / (sc.speed_of_light * sc.epsilon_0)) * au_electric_field_strength_per_SI

    # Create initial state
    rho_0 = np.zeros((n_j, n_j), dtype=complex)
    if temp != 0.0:
        beta = 1/au_energy_per_SI/(sc.Boltzmann*temp)
        g_j = np.diag([n_even if j % 2 == 0 else n_odd for j in range(-j_max, j_max + 1)])
        rho_0 += sl.expm(-beta * B * H_0) @ g_j 
        rho_0 /= np.trace(rho_0)
        if plot:
            plt.imshow(np.abs(rho_0), cmap=cmap)
            plt.show()
    else:
        rho_0[j_max, j_max] = 1.0 + 0.0j
    
    # ODE to solve
    def rho_dot(t, rho):
        H = B*H_0 + H_I * electric_field_squared(t, E0, FWHM) 
        rho = rho.reshape((n_j, n_j))
        rho_dot = -1j * (rho @ H - H @ rho)
        #D_rho = gamma * ( L_op @ rho @ L_op_dag - (L_op_dag @ L_op @ rho + rho @ L_op_dag @ L_op) / 2 )
        #rho_dot += D_rho
        '''for L_op in L_ops:
            L_op_dag = np.conjugate(L_op.T)
            D_rho = gamma * ( L_op @ rho @ L_op_dag - (L_op_dag @ L_op @ rho + rho @ L_op_dag @ L_op) / 2 )
            rho_dot += D_rho'''
        '''
         for i in range(n_j):
            for j in range(n_j):
                rho_dot[i, j] += -(0 if np.abs(i - j) != 2 else pd) * rho[i, j]
        '''

        return rho_dot.flatten()
        #return (-1j * (rho @ H - H @ rho)).flatten()
        

    return solve_ivp(rho_dot, [-3*FWHM, t_end*au_per_ps], rho_0.flatten(), t_eval=np.linspace(-2*FWHM, t_end*au_per_ps, 2000), max_step=5*au_per_ps, method='RK45', rtol=1e-9, atol=1e-9)


# %% -----  SOLVE SIMPLE LvN EQUATION ----- %% #

t1 = time.time()
res = solve_rho(1000, 0.37, 16e9, B, plot=True)
t2 = time.time()
print(f'Time: {t2 - t1:.20f}')
cos2_trace = get_cos2(res)

plt.plot(res.t / au_per_ps, cos2_trace)
plt.figure(figsize=(10, 3))
plt.show()


rho_final = np.reshape(res.y.T[-1], (n_j, n_j))
plt.imshow(np.abs(rho_final), cmap=cmap)

# %% ----- FOCAL VOLUME AVERAGE ----- %% # 
w_pump = 130./2.
w_probe = 26./2.
r_max = 1.7*w_probe 
n_focal = 11

gaussian = lambda r, I0, w: I0 * np.exp(-2 * r**2 / w**2)

# Calculate alignment traces as different spots in the focal volume
def focal_average(I0, B, temp, n_focal, w_pump, w_probe, r_max=w_probe):
    r_lst = np.linspace(0, r_max, n_focal)
    I_pump = gaussian(r_lst, I0, w_pump)
    cos2_lst = []
    for i, (r, I) in enumerate(zip(r_lst, I_pump)):
        cos2_lst.append(get_cos2(solve_rho(1250, temp, I, B)))
        plt.plot(cos2_lst[i])
    cos2_lst = np.array(cos2_lst)
    plt.show()

    # Calculate weights and find total
    weights = r_lst * gaussian(r_lst, 1, w_probe)
    weights /= np.sum(weights)
    cos2_tot = []
    for i in range(len(cos2_lst[0])):
        expval = 0
        for j in range(n_focal):
            expval += weights[j] * cos2_lst[j][i]
        cos2_tot.append(expval)

    return np.array(cos2_tot)

focal_trace = focal_average(15e9, B, 0.37, n_focal, w_pump, w_probe)
# %%

plt.figure(figsize=(6, 4))
plt.ylim(0.45, 0.9)
plt.plot(np.linspace(-2*FWHM / au_per_ps, 1250, len(focal_trace)), focal_trace, 'k-')
ts_pulse = np.linspace(-250, 250, 1000)
plt.fill_between(ts_pulse, electric_field_squared(ts_pulse*au_per_ps, 1, FWHM)*0.2 + 0.49, color='gray', alpha=0.3)
plt.plot(ts_pulse, electric_field_squared(ts_pulse*au_per_ps, 1, FWHM)*0.2 + 0.49, c='gray')
plt.ylabel(r'$\langle\cos^2\theta_\text{2D}\rangle$')
plt.text(500, 0.85, r'$B_\text{2D}=5.25\,\text{GHz}$', fontsize=16)
plt.text(500, 0.82, r'$I_0=23.5\times 10^9$ W/cm$^2$', fontsize=16)
plt.text(500, 0.79, r'$\tau_\text{FWHM}=191$ ps', fontsize=16)
plt.xlabel(r'Delay (ps)')
plt.minorticks_on()
plt.tight_layout()
plt.text(-50, 0.55, r'$I(t)$', fontsize=18, color='grey')
#plt.savefig('Li_adiabatic.png', dpi=350)


# %% ----- PLOT REAL DATA ----- %% #
plt.figure(figsize=(6, 4))
plt.ylim(0.45, 0.9)
ts_pulse = np.linspace(-250, 250, 1000)
xdata, ydata = np.loadtxt('Data/trace1.csv', delimiter=',')
plt.plot(xdata, ydata, 'k.-', lw=1)
offset = 35
plt.text(500, 0.82, r'$I_0=23.5\times 10^9$ W/cm$^2$', fontsize=16)
plt.fill_between(ts_pulse + offset, field(ts_pulse / 0.87)*0.2 + 0.445, color='gray', alpha=0.3)
plt.plot(ts_pulse + offset, field(ts_pulse / 0.87)*0.2 + 0.445, c='gray')
plt.ylabel(r'$\langle\cos^2\theta_\text{2D}\rangle$')
plt.xlabel(r'Delay (ps)')
plt.minorticks_on()
plt.text(-60, 0.5, r'$I(t)$', fontsize=18, color='grey')
plt.tight_layout()
plt.savefig('Li_data1.png', dpi=350)

# %% ----- PLOT THE REST OF THE TRACES ----- %% #
plt.figure(figsize=(6, 4))
plt.ylim(0.47, 0.9)
ts_pulse = np.linspace(-250, 250, 1000)
intensities = [7.2, 10.1, 14.5, 19.5, 23.1]
axsins = plt.gca().inset_axes([0.5, 0.35, 0.45, 0.57])
colors = plt.cm.viridis(np.linspace(0,1,7))
for i, dat in enumerate([4, 6, 2, 3, 5]):
    xdata, ydata = np.loadtxt(f'Data/trace{dat}.csv', delimiter=',')
    plt.plot(xdata - 35, ydata, '-', lw=1, c=colors[i], ms=1, label=intensities[i])
    axsins.plot(xdata - 35, ydata, '.-', lw=1, c=colors[i])
offset = 40
plt.fill_between(ts_pulse + offset, field(ts_pulse)*0.2 + 0.465, color='gray', alpha=0.3)
plt.plot(ts_pulse + offset, field(ts_pulse)*0.2 + 0.465, c='gray')
plt.ylabel(r'$\langle\cos^2\theta_\text{2D}\rangle$')
plt.xlabel(r'Delay (ps)')
plt.minorticks_on()
plt.text(-60, 0.5, r'$I(t)$', fontsize=18, color='grey')
axsins.set(xlim=(200, 900), ylim=(0.48, 0.54))
axsins.minorticks_on()
#axsins.axhline(0.5, c='gray', ls='--')
plt.legend(frameon=False, loc=2, title=r'$I_0$ (GW/cm$^2$)', title_fontsize=14)
plt.xlim(-500, 1100)
plt.tight_layout()
plt.savefig('Li_data2.png', dpi=350)


# %% ----- LINEWIDTH AVERAGE ----- %% #
width_gaussain = lambda B, Delta_B, B_mu: np.exp(-4 * np.log(2) * (B - B_mu)**2 / Delta_B**2)

def linewidth_average(I0, B, temp, n_linewidth, Delta_B):
    Delta_B *= 2*np.pi * 1e6 * au_frequency_per_SI
    #n_linewidth = 11 # ODD!!

    # Calculate the B-consts and their traces
    B_lst = np.linspace(B - Delta_B, B + Delta_B, n_linewidth if n_linewidth % 2 != 0 else n_linewidth + 1)
    cos2_lst = []
    for i, B_i in enumerate(B_lst):
        cos2_lst.append(get_cos2(solve_rho(1250, 0.38, I0, B_i)))
        plt.plot(cos2_lst[i])
    plt.show()

    # Calculate the weights
    weights = width_gaussain(B_lst, Delta_B, B)
    weights /= np.sum(weights)
    cos2_tot = []
    for i in range(len(cos2_lst[0])):
            expval = 0
            for j in range(n_linewidth):
                expval += weights[j] * cos2_lst[j][i]
            cos2_tot.append(expval)
    
    return np.array(cos2_tot)

# %% ----- LINEWIDTH AND INTENSITY ----- %% # 

def linewidth_focal_average(I0, B, temp, n_linewidth, Delta_B, n_focal, w_pump, w_probe, r_max=w_probe):
    Delta_B *= 2*np.pi * 1e6 * au_frequency_per_SI

    # Calculate the B-consts and their traces
    B_lst = np.linspace(B - Delta_B, B + Delta_B, n_linewidth if n_linewidth % 2 != 0 else n_linewidth + 1)
    cos2_lst = []
    for i, B_i in enumerate(B_lst):
        cos2_lst.append(focal_average(I0, B_i, temp, n_focal, w_pump, w_probe, r_max=r_max))
        plt.plot(cos2_lst[i])
    plt.show()

    # Calculate the weights
    weights = width_gaussain(B_lst, Delta_B, B)
    weights /= np.sum(weights)
    cos2_tot = []
    for i in range(len(cos2_lst[0])):
            expval = 0
            for j in range(n_linewidth):
                expval += weights[j] * cos2_lst[j][i]
            cos2_tot.append(expval)
    
    return np.array(cos2_tot)

plt.figure(figsize=(5,2))

cos2_avg = linewidth_focal_average(15e9, 1.29, 0.37, 11, 250, 7, 26, 70)

# %%

plt.figure(figsize=(5.5, 3.5))
times = np.linspace(-2*FWHM/au_per_ps, 1250, 2000)
plt.plot(times, cos2_avg)
plt.ylabel(r'$\langle\cos^2\theta_\text{2D}\rangle$')
plt.xlabel('Time (ps)')
plt.minorticks_on()
plt.tight_layout()

#plt.savefig('Na_simulation.png', dpi=350)


# %%

plt.figure(figsize=(5,2))
plt.plot(linewidth_average(17.5e9, B, 0.37, 11, 250))
plt.axhline(0.5, zorder=-1, c='k')

# %% ----- DROPLET AVERAGE ----- %% #

def droplet_average(I0, B, temp, n_droplet):
    # Calculate angles and their alignent traces
    theta_lst = np.linspace(np.pi/2, 0, n_droplet + 1)[:-1]
    cos2_lst = []
    for theta in theta_lst:
        cos2_lst.append(get_cos2(solve_rho(1250, temp, I0 * np.sin(theta)**2, B)))

    weights = np.abs(np.sin(theta_lst))
    for i, w in enumerate(weights):
        if w != weights[-1]:
            weights[i] *= 2
    weights /= np.sum(weights)
    cos2_tot = [] 
    for i in range(len(cos2_lst[0])):
        expval = 0
        for j in range(n_droplet):
            expval += weights[j] * cos2_lst[j][i]
        cos2_tot.append(expval)
    return cos2_tot

plt.plot(droplet_average(17.5e9, B, 0.37, 20))


# %% -----  DROPLET AND LINEWIDTH AVERAGE ----- %% #

def droplet_linewidth_average(I0, B, temp, n_droplet, n_linewidth, Delta_B):
    theta_lst = np.linspace(np.pi/2, 0, n_droplet + 1)[:-1]
    cos2_lst = []
    for theta in theta_lst:
        cos2_lst.append(linewidth_average(I0*np.sin(theta)**2, B, temp, n_linewidth, Delta_B))

    weights = np.abs(np.sin(theta_lst))
    for i, w in enumerate(weights):
        if w != weights[-1]:
            weights[i] *= 2
    weights /= np.sum(weights)
    cos2_tot = [] 
    for i in range(len(cos2_lst[0])):
        expval = 0
        for j in range(n_droplet):
            expval += weights[j] * cos2_lst[j][i]
        cos2_tot.append(expval)
    return cos2_tot

plt.plot(np.linspace(-2*FWHM/au_per_ps, 1250, 2000), droplet_linewidth_average(17.5e9, B, 0.37, 11, 25, 250))


# %% ----- SOLVE THE SCHRÖDINGER EQUATION ----- %% #
import time

t_start = time.time()
res = solve_ode(800*au_per_ps)
cos2_trace = get_cos2_exp(res)
t_stop = time.time()
print(t_stop - t_start)
plt.plot(res.t / au_per_ps, cos2_trace)


# %% ----- SOLVE THE LvN EQUATION ----- %% #
# Create ODE for LvN equation
def solve_LvN_eqn(T, t_end, fwhm, ratio=0, delta=0):

    # Create inital state
    rho_0 = np.zeros((n_j, n_j), dtype=complex)
    if T == 0:
        rho_0[0, 0] = 1
    else:
        beta = 1/au_energy_per_SI/(sc.Boltzmann*T)
        rho_0 += sl.expm(-beta*H_0)
        rho_0 /= np.trace(rho_0)
        #plt.imshow(np.abs(rho_0))

    def rho_dot(t, rho):
        rho = rho.reshape((n_j, n_j))
        H = H_0 + H_I * (electric_field_squared(t, E0, fwhm))
        return (-1j * (H @ rho - rho @ H)).flatten()
    
    res = solve_ivp(rho_dot, [-2*fwhm, t_end], rho_0.flatten(), max_step=5*au_per_fs, method='RK45', rtol=1e-5, atol=1e-5)
    return res

T = 0.0
res = solve_LvN_eqn(T, 400*au_per_ps, 185*au_per_ps)
cos2_trace = get_cos2_exp(res)
plt.plot(res.t / au_per_ps, cos2_trace)


# %% ----- PLOT THE DENSITY MATRIX ----- %% #

plt.imshow(np.abs(np.reshape(res.y.T[110], (n_j, n_j))), cmap=cmap)
#plt.xlim(20, 80)
#plt.ylim(20, 80)
np.trace(np.reshape(res.y.T[-1], (n_j, n_j)))
#plt.show()
#plt.imshow(H_I)

# %% ----- PLOT POPULATIONS ----- %% #
for i in range(10):
    plt.plot(res.t / au_per_ps, [np.reshape(rho, (n_j, n_j))[i, i] for rho in res.y.T])
plt.xlim(-1, 1)

# %%


n = 1000  # Size of the square 5-diagonal matrix
main_diag = np.random.rand(n)
upper_diag2 = np.random.rand(n-2)
lower_diag2 = np.random.rand(n-2)
B = np.random.rand(n, n)

diags_data = [lower_diag2, np.zeros(n-1), main_diag, np.zeros(n-1), upper_diag2]
offsets = [-2, -1, 0, 1, 2]
A_sparse = sp.diags(diags_data, offsets, shape=(n, n), format='csr')


A_dense = np.diag(main_diag) + np.diag(upper_diag2, 2) + np.diag(lower_diag2, -2)

start_time = time.time()
C = 0.5 * A_sparse @ B
end_time = time.time()
print(f'Sparse: {end_time - start_time:.20f}')

start_time = time.time()
D = 0.5 * A_dense @ B
end_time = time.time()
print(f'Dense: {end_time - start_time:.20f}')

# %%
