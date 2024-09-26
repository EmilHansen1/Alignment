#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "planar_rotor.h"
#include "utils.h"



/**
 * @brief Energy levels of a 2D rotor: Bj²
 * 
 * @param j the angular momentum quantum number
 * @param B the rotational constant
 * @return double 
 */
double E_rot_2D(int j, double B)
{
    return (double) B * j*j;
}


/**
 * @brief Calculates the 2D rotor energy levels
 * The energies lie sequentially as -j, -j+1, ..., j-1, j
 * 
 * @param j_max 
 * @param B 
 * @param rot_energies 
 */
void get_rot_energies_2D(size_t j_max, double B, double rot_energies[j_max + 1])
{
    size_t dim = j_max + 1;
    for(size_t j = 0; j < dim; j++)
    {
        rot_energies[j] = E_rot_2D(j, B);
    }
}


/**
 * @brief Returns the expectation value of cos²θ (in 2D!)
 * 
 * @param j1 first angular momentum quantum number
 * @param j2 second angular monentum quantum number
 * @return double 
 */
double cos2_2D_matelem(int j1, int j2)
{
    if(j1 == j2)
    {
        return 1./2.;
    }
    else if(abs(j1 - j2) == 2)
    {
        return 1./4.;
    }
    else 
    {
        return 0;
    }
}


/**
 * @brief Get the diagonal of the field-free propagator in free 2D rotor basis
 * 
 * @param j_max 
 * @param B 
 * @param dt 
 * @param propagator_diag 
 */
void get_field_free_propagator(int j_max, double B, double dt, dcmplx propagator_diag[j_max+1])
{
    size_t dim = j_max + 1;
    for(size_t j = 0; j < dim; j++)
    {
        propagator_diag[j] = cexp(-I * E_rot_2D(j, B) * dt);
    }
}


/**
 * @brief Performs the field-free propagation
 * 
 * @param j_max
 * @param dt
 * @param n_steps
 * @param time
 * @param B
 * @param psi0
 * @param cos2_exp
*/
void field_free_propagation(const size_t j_max, const double dt, const size_t n_steps, double *time, 
                            const double B, dcmplx psi0[j_max+1], double cos2_exp[n_steps])
{
    // Get the propagator and calculate - only once! - the exponentials
    // Calculation of the complex exponentials becomes extremely time comsuming
    size_t dim = j_max + 1;
    dcmplx U[dim];
    get_field_free_propagator(j_max, B, dt, U);
    dcmplx cos2_psi[dim];   // Placeholder for psi multiplied by cos²

    // Take the rest of the steps
    double current_time = time[0];
    for(size_t i = 0; i < n_steps; i++)
    {
        // Multiply by propagator and update time
        for(size_t j = 0; j < dim; j++)
            psi0[j] *= U[j];
        current_time += dt;
        time[i] = current_time;

        // Calculate alignment expectation value
        //multiply_cos2_2D(j_max, psi0, cos2_psi);
        cos2_exp[i] = get_cos2_2D_expval(j_max, psi0); //(double) scalar_product(dim, psi0, cos2_psi);
    }
}


/**
 * @brief Multiples a vector by the alignemnt cosine matrix in the field-free 2D-rotor representation
 * 
 * @param j_max 
 * @param vec 
 * @param result 
 */
inline void multiply_cos2_2D(const size_t j_max, const dcmplx vec[j_max+1], dcmplx result[j_max+1])
{
    size_t dim = j_max + 1;

    // Set two first entries
    result[0] = vec[0]/2.0 + vec[2]/4.0;
    result[1] = vec[1]/2.0 + vec[3]/4.0;

    // Set the middle dim - 4 entries
    for(size_t i = 2; i < dim - 2; i++)
        result[i] = vec[i-2]/4.0 + vec[i]/2.0 + vec[i+2]/4.0;

    // Set the final two entries - and voíla, we're done!
    result[dim-2] = vec[dim-4]/2.0 + vec[dim-2]/4.0;
    result[dim-1] = vec[dim-3]/2.0 + vec[dim-1]/4.0;
}


double get_cos2_2D_expval(const size_t j_max, const dcmplx psi[j_max+1])
{
    size_t dim = j_max + 1;

    // Keep in mind the degeneracy! It is 2 for all states except j = 0. 
    // We therefore scale and renormalize psi.
    dcmplx psi_scaled[dim];
    for(int i = 0; i < dim; i++)
    {
        double degeneracy = (i == 0) ? 1.0 : 2.0;
        psi_scaled[i] = psi[i] * degeneracy;
    }
    double norm_factor = (double) 1.0/sqrt(scalar_product(dim, psi_scaled, psi_scaled));
    for(int i = 0; i < dim; i++) psi_scaled[i] *= norm_factor;

    // Multiply the state by the matrix representation of cos^2
    dcmplx cos2_psi[dim];
    multiply_cos2_2D(j_max, psi_scaled, cos2_psi);
    return (double) scalar_product(dim, psi_scaled, cos2_psi);
}


/**
 * @brief The TDSE in the format the GSL ode-solver wants
 * 
 * @param t the current time
 * @param _psi wavefunction  at the current time
 * @param _psi_deriv the derivative of the wavefunction at the current time
 * @param params to be cast into ode_params struct
 * @return int 
 */
int planar_rotor_ode(double t, const double _psi[], double _psi_deriv[], void *params)
{
    ode_params *p = (ode_params*) params;
    double field_sq = e_field_squared(t, p->e_field_sq, p->fwhm);

    // GSL requires double[], so we convert back to dcmplx[]
    // Note that _psi_deriv and psi_deriv are just two representations of the SAME array! As for _psi and psi
    const dcmplx *psi = (const dcmplx*) _psi;
    dcmplx *psi_deriv = (dcmplx*) _psi_deriv;

    // Multiply psi by cos² and save in psi_deriv
    multiply_cos2_2D(p->j_max, psi, psi_deriv);

    // The rotor TDSE
    for(int i = 0; i < p->j_max + 1; i++)
        psi_deriv[i] = -I*(p->E_rot[i]*psi[i] - field_sq*p->delta_alpha*psi_deriv[i]/4.0);

    return GSL_SUCCESS;
}


/**
 * @brief Solves the TDSE given the initial state during the interaction with the laser
 * 
 * @param j_max 
 * @param n_steps 
 * @param dt 
 * @param B 
 * @param fwhm 
 * @param e_field_sq 
 * @param delta_alpha 
 * @param time 
 * @param E_rot 
 * @param psi0 
 * @param cos2_exp 
 */
void field_propagation(const size_t j_max, const size_t n_steps, const double dt, const double B, const double fwhm, 
                       const double e_field_sq, const double delta_alpha, double *time, double E_rot[j_max+1], 
                       dcmplx psi0[j_max+1], double cos2_exp[n_steps])
{
    size_t dim = j_max + 1;
    ode_params p = {j_max, e_field_sq, fwhm, E_rot, delta_alpha};
    gsl_odeiv2_system sys = {.dimension = 2*dim, .jacobian=NULL, .function=&planar_rotor_ode, .params=&p};

    // Set up ode driver
    double initial_stepsize = 6.0 * fwhm / ((double) n_steps);
    gsl_odeiv2_driver *ode_driver = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rk8pd, initial_stepsize, 1e-5, 1e-5
    );

    // Variabled for ode-solving
    double ode_timer;       // Time-keeping for the driver
    double t0 = time[0];    // Start time
    dcmplx cos2_psi[dim];   // Placeholder for psi multiplied by cos²

    // Calculate first alignment expectation value
    //multiply_cos2_2D(j_max, psi0, cos2_psi);
    cos2_exp[0] += get_cos2_2D_expval(j_max, psi0); //(double) scalar_product(dim, psi0, cos2_psi);

    // Solve the ode
    ode_timer = t0;
    for(size_t i = 1; i < n_steps; i++)
    {
        gsl_odeiv2_driver_apply(ode_driver, &ode_timer, t0 + (double) i*dt, (double*) psi0);
        time[i] = ode_timer;

        // Calculate alignment expectation value
        //multiply_cos2_2D(j_max, psi0, cos2_psi);
        cos2_exp[i] = get_cos2_2D_expval(j_max, psi0); //(double) scalar_product(dim, psi0, cos2_psi);
        //printf("%f\n", cos2_exp[i]);
    } 
    gsl_odeiv2_driver_free(ode_driver);
}


void propagate_planar_rotor(const size_t j_max, const size_t n_steps_field, const size_t n_steps_field_free )
{

}