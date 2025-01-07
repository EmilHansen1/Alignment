import subprocess


class Simulation:

    def __init__(self, j_max, B, delta_alpha, I0, fwhm, dt, ea, oa, T, t_end, n_lw, dB, w_p, 
                 n_r, n_t, w_x, w_y, custom_flag, custom_fname, homo_abundance, run_type, out_fname):
        
        self.j_max = j_max
        self.B = B
        self.delta_alpha = delta_alpha
        self.I0 = I0
        self.fwhm = fwhm
        self.dt = dt
        self.ea = ea
        self.oa = oa
        self.T = T
        self.t_end = t_end
        self.n_lw = n_lw
        self.dB = dB
        self.w_p = w_p
        self.n_r = n_r
        self.n_t = n_t
        self.w_x = w_x
        self.w_y = w_y
        self.custom_flag = custom_flag
        self.custom_fname = custom_fname
        self.homo_abundance = homo_abundance
        self.run_type = run_type
        self.out_fname = out_fname
        

    def run(self):
        ''' Runs a simulation with the current value for the fields '''

        args = [
            './main',
            str(self.j_max),
            str(self.B),
            str(self.delta_alpha),
            str(self.I0),
            str(self.fwhm),
            str(self.dt),
            str(self.ea),
            str(self.oa),
            str(self.T),
            str(self.t_end),
            str(self.n_lw),
            str(self.dB),
            str(self.w_p),
            str(self.n_r),
            str(self.n_t),
            str(self.w_x),
            str(self.w_y),
            str(self.custom_flag),
            str(self.custom_fname),
            str(self.homo_abundance),
            str(self.run_type),
            str(self.out_fname)
        ]

        subprocess.run(args)



if __name__ == '__main__':
        K2_sim = Simulation(
            j_max=70,
            B=0.69,
            delta_alpha=68.54,
            I0=65.0e9,
            fwhm=0.75,
            dt=0.1,
            ea=5,
            oa=3,
            T=0.37,
            t_end=1250.0,
            n_lw=50,
            dB=50.0,
            w_p=25.0,
            n_r=1,
            n_t=1,
            w_x=70.0,
            w_y=70.0,
            custom_flag=0,
            custom_fname='Pulses/K_pulse_smooth_taper.dat',
            homo_abundance=1.0,
            run_type=0,
            out_fname='cos2.out'
        )

        K2_sim.run()
        