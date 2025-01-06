#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "linear_rotor.h"
#include "utils.h"


/// @brief Function for filling a list with free quantum numbers corresponding to the index
/// @param j_max 
/// @param quantum_numbers A list if free quantum number strutcts
void create_idx_map(const size_t j_max, const int m, free_3d_qn quantum_numbers[j_max])
{
    size_t idx = 0;
    for(size_t j = 0; j < j_max; j++)
    {
        for(int m = -j; m <= j; m++)
        {
            free_3d_qn qn = {j, m};
            quantum_numbers[idx] = qn;
            idx++;
        }
    }
}


/// @brief Helper function for calculating the matrix elements of cos²Θ
/// @param qn The quantum numbers
/// @return The value of a part of the matrix element
double alpha(double j, double m)
{
    return sqrt((j - m + 1.0) * (j + m + 1.0) / ( (2.0 * j + 1.0) * (2.0 * j + 3.0) ));
}


/// @brief Helper function for calculating the matrix elements of cos²Θ
/// @param qn The quantum numbers
/// @return The value of a part of the matrix element
double beta(double j, double m)
{
    return sqrt((j - m) * (j + m) / ( (2.0 * j - 1.0) * (2.0 * j + 1.0) ));
}


/// @brief Calculate the matrix element <jp,mp|cos²Θ|j,m>
/// @param qnp The quantum numbers of the bra
/// @param qn The quantum numbers of the ket
/// @return The value of <jp,mp|cos²Θ|j,m>
double cos2_matelem(free_3d_qn qnp, free_3d_qn qn)
{  
    if(qnp.m == qn.m)
    {
        if (qnp.j == qn.j + 2)
        {
            return alpha(qn.j, qn.m) * alpha(qn.j + 1, qn.m);
        }
        else if(qnp.j == qn.j)
        {
            return alpha(qn.j, qn.m) * beta(qn.j + 1, qn.m) + alpha(qn.j - 1, qn.m) * beta(qn.j, qn.m);
        }
        else if(qnp.j == qn.j - 2)
        {
            return beta(qn.j, qn.m) * beta(qn.j - 1, qn.m);
        }
    }
    return 0.0;
}



/// @brief Fills up the diagonals of the tridiag struct
/// @param j_max The maximum j quantum number (not included)
/// @param quantum_numbers The index -> quantum numbers map (array)
/// @param diags The tridiags struct to be filled
void fill_diagonals(const size_t j_max, const free_3d_qn quantum_numbers[j_max * j_max], tridiag *diags)
{
    const size_t dim = j_max * j_max;

    // First fill up the diagnonal
    for(size_t i = 0; i < dim; i++)
    {
        diags->diag[i] = cos2_matelem(quantum_numbers[i], quantum_numbers[i]);
    }

    // Then the off-diagonal
    for(size_t i = 0; i < dim - 2; i++)
    {
        diags->off_diag[i] = cos2_matelem(quantum_numbers[i], quantum_numbers[i+2]);
    }
}


/// @brief Multiply the tri-diagonal cos²Θ matrix on the wavefunction, given its diagonals
/// @param j_max The maximum j quantum number (not included)
/// @param cos2_theta_psi Vector to hold the result
/// @param psi The wavefunction
/// @param diag The diagonal of the matrix representation (Δj = 0)
/// @param offdiag The off-diagonal of the matrix representation (Δj = ±2)
void multiply_cos2_3D(const size_t j_max, dcmplx cos2_theta_psi[j_max * j_max], const dcmplx psi[j_max * j_max], 
                      const double diag[j_max * j_max], const double offdiag[j_max * j_max - 2])
{
    const size_t dim = j_max * j_max;

    // We start with the first two entries
    cos2_theta_psi[0] = diag[0] * psi[0] + offdiag[0] * psi[2];
    cos2_theta_psi[1] = diag[1] * psi[1] + offdiag[1] * psi[3];

    // Then loop over all except the last (and first) two indices
    for(size_t i = 2; i < dim - 2; i++)
        cos2_theta_psi[i] = diag[i] * psi[i] + offdiag[i - 2] * psi[i - 2] + offdiag[i] * psi[i + 2];

    // Finally the last two
    cos2_theta_psi[dim - 2] = diag[dim - 2] * psi[dim - 2] + offdiag[dim - 4] * psi[dim - 4];
    cos2_theta_psi[dim - 1] = diag[dim - 1] * psi[dim - 1] + offdiag[dim - 3] * psi[dim - 3];
}




/// @brief The energy of a free 3D rotor given the quantum numbers and the rotational constant
/// @param qn The qunatum number struct
/// @param B The rotational constant
/// @return The energy of a free 3D rotor
double E_rot_3D(const free_3d_qn qn, const double B)
{
    return B * qn.j * (qn.j + 1.0);
}


/// @brief Fills an array with the free rotor Hamiltonian eigenvalues (with degeneracy)
/// @param j_max The maximum quantum number (not included)
/// @param B The rotational constant
/// @param energies The list of energies to be filled
void get_rot_energies_3D(const size_t j_max, const double B, double energies[j_max * j_max], const free_3d_qn quantum_numbers[j_max * j_max])
{
    const size_t dim = j_max * j_max;
    for(size_t i = 0; i < dim; i++)
    {
        energies[i] = E_rot_3D(quantum_numbers[i], B);
    }
}


/// @brief Produces the field-free propagator for the post-pulse dynamics
/// @param j_max The maximum j quantum number (not included)
/// @param B The rotational constant
/// @param dt The timestep
/// @param quantum_numbers The index to quantum number map 
/// @param propagator_diag The diagonal of the propagator to be filled
void get_field_free_prpagator_3D(const size_t j_max, const double B, double dt, const free_3d_qn quantum_numbers[j_max * j_max], dcmplx propagator_diag[j_max * j_max])
{
    const size_t dim = j_max * j_max;
    double energies[dim];
    get_rot_energies_3D(j_max, B, energies, quantum_numbers);
    for(size_t i = 0; i < dim; i++)
    {
        propagator_diag[i] = cexp(-I * energies[i] * dt); 
    }
}


/// @brief 
/// @param j_max 
/// @param quantum_numbers 
/// @param dt 
/// @param n_steps 
/// @param counter 
/// @param B 
/// @param psi0 
/// @param cos2_exp 
/// @param weight 
/// @param diags 
void field_free_propagation_3D(const size_t j_max, const free_3d_qn quantum_numbers[j_max * j_max], const double dt, const size_t n_steps, 
                               size_t *counter, const double B, dcmplx psi0[j_max * j_max], double cos2_exp[], double weight, const tridiag *diags)
{
    // Get field-free propagator
    const size_t dim = j_max * j_max;
    dcmplx U_diag[dim];
    get_field_free_prpagator_3D(j_max, B, dt, quantum_numbers, U_diag);
    dcmplx cos2_psi[dim];

    // Take the steps
    for(size_t i = 0; i < n_steps; i++)
    {
        // Propagate the wavefunction
        for(size_t j = 0; j < dim; j++)
        {
            psi0[i] *= U_diag[i];
        }

        // Calculate degree of alignment
        multiply_cos2_3D(j_max, cos2_psi, psi0, diags->diag, diags->off_diag);
        cos2_exp[*counter] += (double) weight * scalar_product(dim, psi0, cos2_psi);
        (*counter)++;
    }
}



/// @brief 
/// @param t 
/// @param _psi 
/// @param _psi_deriv 
/// @param params 
/// @return 
int linear_rotor_ode(double t, const double _psi[], double _psi_deriv[], void *params)
{
    const ode_params *p = (ode_params*) params;
    const int dim = p->j_max * p->j_max;
    const double field_sq = (p->laser->custom_pulse_flag) 
                            ? e_field_squared_custom(t, p->e_field_sq, p->laser->spline, p->laser->acc) 
                            : e_field_squared(t, p->e_field_sq, p->fwhm);

    // GSL requires double[], so we convert back to dcmplx[]
    // Note that _psi_deriv and psi_deriv are just two representations of the SAME array! As for _psi and psi
    const dcmplx *psi = (const dcmplx*) _psi;
    dcmplx *psi_deriv = (dcmplx*) _psi_deriv;
    
    // Multiply psi by cos² and save in psi_deriv
    multiply_cos2_3D(p->j_max, psi_deriv, psi, p->diagonals->diag, p->diagonals->off_diag);

    // The rotor TDSE
    for(size_t i = 0; i < dim; i++)
        psi_deriv[i] = -I*(p->E_rot[i]*psi[i] - field_sq*p->delta_alpha*psi_deriv[i]/4.0);

    return GSL_SUCCESS;
}



void field_propagation_3D(const size_t j_max, const size_t n_steps, const double dt, const double B, const double fwhm, 
                          const double e_field_sq, const double delta_alpha, const double E_rot[2*j_max+1], 
                          dcmplx psi0[2*j_max+1], double t0, size_t *counter, double *cos2_exp, double weight, const field_params *fp, const tridiag *diags)
{
    const size_t dim = j_max * j_max;
    ode_params p = {j_max, e_field_sq, fwhm, E_rot, delta_alpha, fp, diags};
    gsl_odeiv2_system sys = {.dimension=2*dim, .jacobian=NULL, .function=&linear_rotor_ode, .params=&p};

    // Set up ode driver
    double initial_stepsize = 40.0 * (1e-12 * 4.134137333518211e+16) / ((double) n_steps);
    gsl_odeiv2_driver *ode_driver = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rk8pd, initial_stepsize, 1e-5, 1e-5
    );

    // Variabled for ode-solving
    double ode_timer;       // Time-keeping for the driver
    dcmplx cos2_psi[dim];   // Placeholder for psi multiplied by cos²

    // Calculate first alignment expectation value
    multiply_cos2_3D(j_max, psi0, cos2_psi, p.diagonals->diag, p.diagonals->off_diag);
    cos2_exp[*counter] += (double) weight * scalar_product(dim, psi0, cos2_psi);
    (*counter)++;

    // Solve the ode
    ode_timer = t0;
    for(size_t i = 1; i < n_steps; i++)
    {
        // Apply the driver
        gsl_odeiv2_driver_apply(ode_driver, &ode_timer, t0 + (double) i*dt, (double*) psi0);

        // Calculate alignment expectation value
        multiply_cos2_3D(j_max, psi0, cos2_psi, p.diagonals->diag, p.diagonals->off_diag);
        cos2_exp[*counter] += (double) weight * scalar_product(dim, psi0, cos2_psi);
        (*counter)++;
    } 
    gsl_odeiv2_driver_free(ode_driver);
}

void linear_rotor_propagation(solver_params *params, dcmplx psi0[params->dim], double cos2[], double weight)
{
    // Counter for keeping track of the time(step)
    size_t counter = 0;

    // Solve the ODE during the field
    field_propagation_3D(params->molecule->j_max, 
                         params->n_steps_field,
                         params->dt,
                         params->molecule->B,
                         params->field->fwhm,
                         params->field->e_field_squared,
                         params->molecule->delta_alpha,
                         params->molecule->E_rot,
                         psi0,
                         params->t_start,
                         &counter,
                         cos2,
                         weight, 
                         params->field,
                         params->diagonals);

    // Solve the field-free evolution 
    /*field_free_propagation_3D(params->molecule->j_max,
                              params->molecule->qn, 
                              params->dt,
                              params->n_steps_field_free,
                              &counter,
                              params->molecule->B,
                              psi0,
                              cos2,
                              weight,
                              params->diagonals);*/


}