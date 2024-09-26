#ifndef PLANAR_ROTOR_H
#define PLANAR_ROTOR_H

typedef double _Complex dcmplx;

typedef struct 
{
    const size_t j_max;         // The maximum angular momentum quantum number
    const double e_field_sq;    // The peak E-field squared
    const double fwhm;          // The full width at half maximum of the of the Gaussian pulse 
    const double *E_rot;        // Array with the rotational energy levels (from -j to j)
    const double delta_alpha;   // Polarizability anisotropy
} ode_params; 

double E_rot_2D(int j, double B);

void get_rot_energies_2D(size_t j_max, double B, double rot_energies[j_max+1]);

void get_field_free_propagator(int j_max, double B, double dt, dcmplx propagator_diag[j_max+1]);

void multiply_cos2_2D(const size_t j_max, const dcmplx vec[j_max+1], dcmplx result[j_max+1]);

double get_cos2_2D_expval(const size_t j_max, const dcmplx psi[j_max+1]);

void field_free_propagation(const size_t j_max, const double dt, const size_t n_steps, double *time, 
                            const double B, dcmplx psi0[j_max+1], double cos2_exp[n_steps]);

int planar_rotor_ode(double t, const double _psi[], double _psi_deriv[], void* params);

void field_propagation(const size_t j_max, const size_t n_steps, const double dt, const double B, const double fwhm, 
                       const double e_field_sq, const double delta_alpha, double *time, double E_rot[j_max+1], 
                       dcmplx psi0[j_max+1], double cos2_exp[n_steps]);

#endif  /* PLANAR_ROTOR_H */