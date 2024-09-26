# %%
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.sparse as ss
from scipy.sparse.linalg import expm, expm_multiply
import scipy.linalg as sl
from scipy.integrate import solve_ivp
import time
from multiprocessing import Pool

# Matplotlib asthetics
major = 6 # 6
minor = 3 # 3
width = 1.5
plt.rc('text', usetex=False)
#plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
plt.rc("axes", labelsize=16) # 18
plt.rc("xtick", labelsize=14, top=True, direction="in")
plt.rc("ytick", labelsize=14, right=True, direction="in")
plt.rc("axes", titlesize=18)
plt.rc("legend", fontsize=14)
plt.rcParams['font.family'] = 'sans-serif'
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

# %%


class PlanarRotor:

    def __init__(self, B: float, delta_alpha: float, T: float, j_max: int, eo_pair: tuple[int, int]=(1, 1)):
        # Fields
        self._B = B
        self.delta_alpha = delta_alpha
        self.T = T
        self.j_max = j_max
        self.eo_pair = eo_pair
        self.n = 2 * self.j_max + 1

        # Solution objects
        self.solution = None
        self.trace = None
        self.times = None

        # Allocate matrix representations
        self.cos2_theta_mat = np.zeros((self.n, self.n))
        self.J_squared_mat = np.zeros((self.n, self.n))
        self.cos2_theta_mat_s = np.zeros((self.n, self.n))
        self.J_squared_mat_s = np.zeros((self.n, self.n))
        self.calculate_operators()
        

    @property
    def B(self):
        return self._B
    

    @property
    def B_au(self):
        return self._B * (2. * np.pi) * 1.e9 *au_frequency_per_SI
    

    @property
    def delta_alpha_au(self):
        return self.delta_alpha * au_per_angstrom**3
    

    @B.setter
    def B(self, new_B):
        self._B = new_B


    def energy(self, j):
        assert abs(j) < self.j_max - 1, 'j larger than j_max-1'
        return self.B * j**2
    

    def envelope_squared(self, t, amplitude, fwhm):
        return amplitude**2 * np.exp(-four_log_two * (t / fwhm)**2)
        

    def _gaussian_peak_intensity(self, E, tau, wx, wy):
        return 4 * np.sqrt(np.log(2)) * E / (np.pi**(3/2) * tau * wx * wy)
    

    def _intensity_to_amplitude(self, intensity_si):
        amplitude = np.sqrt(2 * (intensity_si * 1e4) / (sc.speed_of_light * sc.epsilon_0)) * au_electric_field_strength_per_SI
        return amplitude


    def _cos2_squared_matelem(self, j1, j2):
        if np.abs(j1 - j2) == 2:
            return 1./4.
        elif np.abs(j1 - j2) == 0:
            return 1./2.
        else:
            return 0
    

    def calculate_operators(self):
        # First the non-sparse matrices
        j_lst = np.arange(self.n) - self.j_max
        for i, j1 in enumerate(j_lst):
            self.J_squared_mat[i,i] = j1**2
            for j, j2 in enumerate(j_lst):
                self.cos2_theta_mat[i,j] = self._cos2_squared_matelem(j1, j2)

        # And then sparsify
        self.J_squared_mat_s = ss.csr_matrix(self.J_squared_mat.copy())
        self.cos2_theta_mat_s = ss.csr_matrix(self.cos2_theta_mat.copy())


    def calculate_cos2_exp(self):
        # Calculate the degree of alignment from a ODE solution object
        cos2_expval = []
        for rho_arr in self.solution.y.T:
            rho = np.reshape(rho_arr, (self.n, self.n))
            cos2_expval.append(np.real(np.trace(self.cos2_theta_mat @ rho)))
        return np.array(cos2_expval)


    def _get_thermal_state(self):
        # Creates a thermal density matrix with the proper nuclear-spin statistics
        rho_0 = np.zeros((self.n, self.n), dtype=complex)
        beta = 1/au_energy_per_SI/(sc.Boltzmann*self.T)
        n_even, n_odd = list(self.eo_pair)
        j_range = range(-self.j_max, self.j_max + 1)
        g_j = np.diag([n_even if j % 2 == 0 else n_odd for j in j_range])
        rho_0 += sl.expm(-beta * self.B_au * self.J_squared_mat) @ g_j 
        rho_0 /= np.trace(rho_0)
        return rho_0
    

    def solve_thermal_ode(self, t_range, peak_intensity, fwhm, plot=False, n_eval=1000):
        # Create the initial state
        rho = self._get_thermal_state()
        amplitude = self._intensity_to_amplitude(peak_intensity)
        envelope_sq = lambda t: self.envelope_squared(t, amplitude, fwhm)

        # Define the ODE
        def _ode(t, rho):
            H = self.B_au * self.J_squared_mat - envelope_sq(t) * self.delta_alpha_au * self.cos2_theta_mat / 4.0
            rho = rho.reshape((self.n, self.n))
            rho_dot = 1.j * (rho @ H - H @ rho)
            return rho_dot.flatten()

        # Solve the ODE
        print('\n===== SOLVING ODE =====')
        t_eval = np.linspace(*t_range, n_eval)
        self.times = t_eval.copy()
        t_start = time.time()
        self.solution = solve_ivp(_ode, t_range, rho.flatten(), t_eval=t_eval, max_step=5*au_per_ps)
        t_end = time.time()
        print(f'ODE solved in: {t_end - t_start:.3f} s\n')

        # Calculate expectation values
        print('\n===== CALCULATING EXPECTATION VALUES =====')
        t_start = time.time()
        self.trace = self.calculate_cos2_exp()
        t_end = time.time()
        print(f'Expectation values calculated in: {t_end - t_start:.3f} s\n')

    
    def get_expectation_value(self, rho, O):
        return np.real(np.trace(rho @ O))

    def solve_thermal_direct(self, t_range, peak_intensity, fwhm, n_eval=1000, dt=None, weight=1., B=None):
        # Create the initial state
        rho = self._get_thermal_state()
        amplitude = self._intensity_to_amplitude(peak_intensity)
        envelope_sq = lambda t: self.envelope_squared(t, amplitude, fwhm)

        if B is not None:
            self._B = B

        # Times
        t0 = t_range[0]
        if dt is None:
            dt = int((t_range[1] - t_range[0]) / n_eval)
        times, expval = [], []

        # Field-free time-evolution operators
        U_ff = ss.dia_matrix(expm(-1.j * self.B_au * self.J_squared_mat_s * dt), dtype=complex)
        U_ff_H = U_ff.H

        # Iterate
        t_start = time.time()
        for i in range(n_eval):
            t = t0 + i * dt
            times.append(t)

            # Calculate expectation value
            expval.append(self.get_expectation_value(rho, self.cos2_theta_mat_s))

            # Propagate a single step
            if np.abs(t) < 3 * fwhm:
                H = self.B_au * self.J_squared_mat_s - envelope_sq(t) * self.delta_alpha_au * self.cos2_theta_mat_s / 4.0
                rho = expm_multiply(-1.j * H * dt, rho) @ expm(1.j * H * dt)
            else: 
                rho = U_ff @ rho @ U_ff_H
        t_end = time.time()

        print(f'DIRECT SOLUTION DONE IN: {t_end - t_start:.3f}')
        print(f'Trace: {np.trace(rho)}')
        self.trace = weight * np.array(expval)
        self.times = np.array(times)
        return self.trace.copy(), self.times.copy()


    def calculate_total_average(self, t_range, peak_intensity, fwhm, n_focal, w_pump, w_probe, n_incoh, Delta_B, plot=False, 
                                n_eval=1000, method='ode'):
        # Intensities and weights for focal averaging
        gaussian = lambda r, I0, w: I0 * np.exp(-2 * r**2 / w**2)
        r_lst = np.linspace(0.01*w_probe, 1.7*w_probe, n_focal)
        I_pump = gaussian(r_lst, peak_intensity, w_pump)
        I_weights = r_lst * gaussian(r_lst, 1, w_probe)
        I_weights /= np.sum(I_weights)

        # B-constants and weights for incoherent averaging
        # Make sure n_incoh is odd
        width_gaussian = lambda B, Delta_B, B_mu: np.exp(-4 * np.log(2) * (B - B_mu)**2 / Delta_B**2)
        Delta_B *= 1e-3
        n_incoh = n_incoh if n_incoh % 2 != 0 else n_incoh + 1
        B_lst = np.linspace(self.B - Delta_B, self.B + Delta_B, n_incoh)
        B_weights = width_gaussian(B_lst, Delta_B, self.B)
        B_weights /= np.sum(B_weights)

        # ... Loop over it all
        traces = np.zeros((n_focal, n_incoh, n_eval))
        norm_factor = 1/( np.sum(I_weights) * np.sum(B_weights))


        '''arguments = [(intensity, B_val, w_I*w_B/norm_factor) 
                     for intensity, B_val, w_I, w_B in zip(I_pump, B_lst, I_weights, B_weights)]
        
        def exec_func(I, B, w):
            return self.solve_thermal_direct(t_range, I, fwhm, n_eval=n_eval, B=B, weight=w)
        
        with Pool(processes=1) as pool:
            trace = pool.starmap(exec_func, arguments)
        return np.array(trace)'''

        for i, (w_I, intensity) in enumerate(zip(I_weights, I_pump)):
            for j, (w_B, rot_const) in enumerate(zip(B_weights, B_lst)):
                # Solve and calculate
                self.B = rot_const
                if method == 'ode':
                    self.solve_thermal_ode(t_range, intensity, fwhm)
                    traces[i, j] = w_I * w_B * self.calculate_cos2_exp() / norm_factor
                else:
                    trace, _ = self.solve_thermal_direct(t_range, intensity, fwhm, n_eval=n_eval)
                    traces[i, j] = w_I * w_B * trace.copy() / norm_factor
        self.traces = traces

# %% Na simulation
B = 1.33    # GHz
Δα = 31.2   # Å³
I0 = 50e9   # GW/cm²
T = 0.37     # K
FWHM = 1 * au_per_ps
t_range = [-2*FWHM, 1250*au_per_ps]
j_max = 100 
eo_pair = (5, 3)

# Solve single ensemble
rotor = PlanarRotor(B, Δα, T, j_max, eo_pair=eo_pair)
#rotor.solve_thermal_ode(t_range, I0, FWHM)
#plt.plot(rotor.times / au_per_ps, rotor.trace)
#plt.show()

# Include all averaging
n_focal, n_incoh = 1, 11
w_pump = 40.
w_probe = 25.
Delta_B = 100 # MHz
trace_d = rotor.solve_thermal_direct(t_range, I0, FWHM, n_eval=2000)
#plt.plot(times_d / au_per_ps, trace_d)


rotor.calculate_total_average(t_range, I0, FWHM, n_focal, w_pump, w_probe, n_incoh, Delta_B, 
                              n_eval=2000, method='direct')    

print(rotor._gaussian_peak_intensity(170e-6, 171e-12, 179.4e-6/2, 80.54e-6/2) / 1e9 * 1e-4)

traces = np.reshape(rotor.traces, (n_focal * n_incoh, len(rotor.times)))
Na_trace = np.sum(traces, axis=0)
plt.plot(rotor.times / au_per_ps, np.sum(traces, axis=0))

# %% Li simulation
B = 5.25    # GHz
Δα = 65.73   # Å³
I0 = 23.5e9   # GW/cm²
T = 0.37     # K
FWHM = 191 * au_per_ps
t_range = [-2*FWHM, 1250*au_per_ps]
j_max = 10 
eo_pair = (5, 3)

# Solve single ensemble
rotor = PlanarRotor(B, Δα, T, j_max, eo_pair=eo_pair)

# Include all averaging
n_focal, n_incoh = 11, 11
w_pump = 40.
w_probe = 26.
Delta_B = 200 # MHz
rotor.calculate_total_average(t_range, I0, FWHM, n_focal, w_pump, w_probe, n_incoh, Delta_B, n_eval=2000) 

print(rotor._gaussian_peak_intensity(170e-6, 171e-12, 179.4e-6/2, 80.54e-6/2) / 1e9 * 1e-4)

traces = np.reshape(rotor.traces, (n_focal * n_incoh, len(rotor.times)))
Li_trace = np.sum(traces, axis=0)
plt.plot(rotor.times / au_per_ps, np.sum(traces, axis=0))

# %% Na simulation
B = 1.33    # GHz
Δα = 31.2   # Å³
I0 = 16.459e9   # GW/cm²
T = 0.37     # K
FWHM = 171 * au_per_ps
t_range = [-2*FWHM, 1250*au_per_ps]
j_max = 12 
eo_pair = (5, 3)

# Solve single ensemble
rotor = PlanarRotor(B, Δα, T, j_max, eo_pair=eo_pair)
rotor.solve_thermal_ode(t_range, I0, FWHM)
print(rotor.B_au, rotor.B)
#plt.plot(rotor.times / au_per_ps, rotor.trace)
#plt.show()

# Include all averaging
n_focal, n_incoh = 11, 11
w_pump = 40.
w_probe = 25.
Delta_B = 200 # MHz
rotor.calculate_total_average(t_range, I0, FWHM, n_focal, w_pump, w_probe, n_incoh, Delta_B, n_eval=2000, method='direct')    

print(rotor._gaussian_peak_intensity(170e-6, 171e-12, 179.4e-6/2, 80.54e-6/2) / 1e9 * 1e-4)

traces = np.reshape(rotor.traces, (n_focal * n_incoh, len(rotor.times)))
Na_trace = np.sum(traces, axis=0)
plt.plot(rotor.times / au_per_ps, np.sum(traces, axis=0))

# %% K simulation
B = 0.7 #0.89    # GHz
Δα = 68.54   # Å³
I0 = 11.78e9   # W/cm²
T = 0.37     # K
FWHM = 185 * au_per_ps
t_range = [-2*FWHM, 1250*au_per_ps]
j_max = 30 
eo_pair = (5, 3)
rotor = PlanarRotor(B, Δα, T, j_max, eo_pair=eo_pair)

# Include all averaging
n_focal, n_incoh = 3, 3
w_pump = 40.
w_probe = 26.
Delta_B = 400 # MHz
rotor.calculate_total_average(t_range, I0, FWHM, n_focal, w_pump, w_probe, n_incoh, Delta_B, 
                              n_eval=700, method='direct')    

traces = np.reshape(rotor.traces, (n_focal * n_incoh, len(rotor.times)))
K_trace = np.sum(traces, axis=0)
plt.plot(rotor.times / au_per_ps, np.sum(traces, axis=0))

# %% Rb simulation
B = 0.28 #0.89    # GHz
Δα = 67.82   # Å³
I0 = 18.7e9    # W/cm² #rotor._gaussian_peak_intensity(160e-6, 179e-12, 63.73e-6/2, 179.0e-6/2) * 1e-4 / 1e9
T = 0.37     # K
FWHM = 179 * au_per_ps
t_range = [-2*FWHM, 1250*au_per_ps]
j_max = 60 
eo_pair = (7, 5)
rotor = PlanarRotor(B, Δα, T, j_max, eo_pair=eo_pair)

# Include all averaging
n_focal, n_incoh = 3, 3
w_pump = 40.
w_probe = 26.
Delta_B = 50 # MHz
rotor.calculate_total_average(t_range, I0, FWHM, n_focal, w_pump, w_probe, n_incoh, Delta_B, n_eval=600, method='direct')

traces = np.reshape(rotor.traces, (n_focal * n_incoh, len(rotor.times)))
Rb_trace = np.sum(traces, axis=0)
plt.plot(rotor.times / au_per_ps, np.sum(traces, axis=0))

# %% Create plot

fig, axs = plt.subplots(4, 1, figsize=(6, 8), sharex='all', constrained_layout=True)#,  gridspec_kw=dict(hspace=0, wspace=0))
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0)
fig.get_layout_engine().set(w_pad=0, h_pad=0, hspace=0, wspace=0)

y_data = [Li_trace, Na_trace, K_trace, Rb_trace]
FWHMs = [191, 171, 185, 179]
I0s = [23.5, 16.4, 11.7, 18.7]
target = [r'Li$_2$', r'Na$_2$', r'K$_2$', r'Rb$_2$']
x_data = rotor.times
colors = ['red', 'blue', 'darkblue', 'purple']
#colors = ['silver', 'darkgrey', 'dimgrey', 'black']
#cmap = mpl.colormaps['RdBu']
#colors = cmap([0.1, 0.3, 0.7, 0.9])

for i, ax in enumerate(axs):
    times = np.linspace(-2*FWHMs[i], 1250, len(y_data[i]))
    ax.plot(times, y_data[i], color=colors[i], lw=2)
    #ax.minorticks_on()
    ax.set_ylabel(r'$\langle\cos^2\theta_\text{2D}\rangle$')
    ax.text(0.97, 0.8, f'{target[i]}$\,(1^3\Sigma_\mathrm{{u}}^+)$', fontsize=14, transform=ax.transAxes, ha='right')
    ax.text(0.97, 0.65, f'$I_0=\\mathrm{{{I0s[i]}}}$ GW/cm$^2$', fontsize=14, transform=ax.transAxes, ha='right')
    ax.text(0.97, 0.5, f'$\\tau={FWHMs[i]}$ $\\mathrm{{ps}}$', fontsize=14, transform=ax.transAxes, ha='right')
    ax.text(0.97, 0.35, f'$\\Delta B={200}$ MHz', fontsize=14, transform=ax.transAxes, ha='right')
axs[-1].set_xlabel(r'$t$ (ps)')
#plt.tight_layout()
plt.savefig('2D_Ak2_rotation.png', dpi=500)