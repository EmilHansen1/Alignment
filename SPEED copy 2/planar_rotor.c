#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "planar_rotor.h"
#include "utils.h"
//#include "units.h"



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
void get_rot_energies_2D(size_t j_max, double B, double rot_energies[2*j_max + 1])
{
    int j = -j_max;
    size_t dim = 2*j_max + 1;
    for (size_t i = 0; i < dim; i++)
    {
        rot_energies[i] = E_rot_2D(j, B);
        j++;
    }
}


/**
 * @brief Returns the expectation value of cos²θ (in 2D!)
 * 
 * @param j1 first angular momentum quantum number
 * @param j2 second angular monentum quantum number
 * @return double <j1|cos²θ|j2>
 */
double cos2_2D_matelem(int j1, int j2)
{
    if (j1 == j2)
    {
        return 1./2.;
    }
    else if (abs(j1 - j2) == 2)
    {
        return 1./4.;
    }
    else 
    {
        return 0;
    }
}


/**
 * @brief Calculate all exponentials for FF-propagation and store in a matrix
 * 
 * @param j_max 
 * @param B 
 * @param dt 
 * @param exponentials 
 */
void get_exponentials(int j_max, double B, double dt, dcmplx exponentials[2*j_max+1][2*j_max+1])
{
    int j1 = -j_max;
    int j2 = -j_max;
    for (size_t i = 0; i < 2*j_max+1; i++)
    {
        for (size_t j = 0; j < 2*j_max+1; j++)
        {
            exponentials[i][j] = cexp(-I * B * ((double)(j2*j2 - j1*j1)) * dt);
            j2++;
        }
        j2 = -j_max;
        j1++;
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
void get_field_free_propagator(int j_max, double B, double dt, dcmplx propagator_diag[2*j_max+1])
{
    int j = -j_max;
    for (size_t i = 0; i < 2*j_max+1; i++)
    {
        propagator_diag[i] = cexp(-I * B * ((dcmplx) j*j) * dt);
        j++;
    }
}


void field_free_propagation(const size_t j_max, const double dt, const size_t n_steps, size_t *counter, 
                            const double B, dcmplx psi0[2*j_max+1], double cos2_exp[], double weight)
{
    // Get the propagator and calculate - only once! - the exponentials
    // Calculation of the complex exponentials becomes extremely time comsuming
    size_t dim = 2*j_max + 1;
    dcmplx U[dim];
    get_field_free_propagator(j_max, B, dt, U);
    dcmplx cos2_psi[dim];   // Placeholder for psi multiplied by cos²

    // Take the rest of the steps
    for (size_t i = 0; i <= n_steps; i++)
    {
        // Multiply by propagator and update time
        for (size_t j = 0; j < dim; j++)
            psi0[j] *= U[j];

        // Calculate alignment expectation value
        multiply_cos2_2D(j_max, psi0, cos2_psi);
        cos2_exp[*counter] += (double) weight * scalar_product(dim, psi0, cos2_psi);
        //printf("%f\n", cos2_exp[*counter]);
        (*counter)++;
    }
}


/**
 * @brief Multiples a vector by the alignemnt cosine matrix in the field-free 2D-rotor representation
 * 
 * @param j_max 
 * @param vec 
 * @param result 
 */
inline void multiply_cos2_2D(const size_t j_max, const dcmplx vec[2*j_max+1], dcmplx result[2*j_max+1])
{
    size_t dim = 2*j_max + 1;

    // Set two first entries
    result[0] = vec[0]/2.0 + vec[2]/4.0;
    result[1] = vec[1]/2.0 + vec[3]/4.0;

    // Set the middle dim - 4 entries
    for (size_t i = 2; i < dim - 2; i++)
        result[i] = vec[i-2]/4.0 + vec[i]/2.0 + vec[i+2]/4.0;

    // Set the final two entries - and voíla, we're done!
    result[dim-2] = vec[dim-4]/2.0 + vec[dim-2]/4.0;
    result[dim-1] = vec[dim-3]/2.0 + vec[dim-1]/4.0;
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
    double field_sq = (p->laser->custom_pulse_flag) ? 
                      e_field_squared_custom(t, p->e_field_sq, p->laser->spline, p->laser->acc) : e_field_squared(t, p->e_field_sq, p->fwhm);

    // GSL requires double[], so we convert back to dcmplx[]
    // Note that _psi_deriv and psi_deriv are just two representations of the SAME array! As for _psi and psi
    const dcmplx *psi = (const dcmplx*) _psi;
    dcmplx *psi_deriv = (dcmplx*) _psi_deriv;
    
    // Multiply psi by cos² and save in psi_deriv
    multiply_cos2_2D(p->j_max, psi, psi_deriv);

    // The rotor TDSE
    for (int i = 0; i < 2*p->j_max + 1; i++)
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
                       const double e_field_sq, const double delta_alpha, const double E_rot[2*j_max+1], 
                       dcmplx psi0[2*j_max+1], double t0, size_t *counter, double *cos2_exp, double weight, const field_params *fp)

{
    size_t dim = 2*j_max + 1;
    ode_params p = {j_max, e_field_sq, fwhm, E_rot, delta_alpha, fp};
    gsl_odeiv2_system sys = {.dimension = 2*dim, .jacobian=NULL, .function=&planar_rotor_ode, .params=&p};

    // Set up ode driver
    double initial_stepsize = 4.0 * fwhm * (1e-12 * 4.134137333518211e+16) / ((double) n_steps);
    gsl_odeiv2_driver *ode_driver = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rk8pd, initial_stepsize, 1e-5, 1e-5
    );

    // Variabled for ode-solving
    double ode_timer;       // Time-keeping for the driver
    dcmplx cos2_psi[dim];   // Placeholder for psi multiplied by cos²

    // Calculate first alignment expectation value
    multiply_cos2_2D(j_max, psi0, cos2_psi);
    cos2_exp[*counter] += (double) weight * scalar_product(dim, psi0, cos2_psi);
    (*counter)++;

    // Solve the ode
    ode_timer = t0;
    for (size_t i = 1; i < n_steps; i++)
    {
        // Apply the driver
        gsl_odeiv2_driver_apply(ode_driver, &ode_timer, t0 + (double) i*dt, (double*) psi0);

        // Calculate alignment expectation value
        multiply_cos2_2D(j_max, psi0, cos2_psi);
        cos2_exp[*counter] += (double) weight * scalar_product(dim, psi0, cos2_psi);
        (*counter)++;
    } 
    gsl_odeiv2_driver_free(ode_driver);
}


void planar_rotor_propagation(solver_params *params, dcmplx psi0[params->dim], double cos2[], double weight)
{
    // Counter for keeping track of the time(step)
    size_t counter = 0;

    // Solve the ODE during the field
    field_propagation(params->molecule.j_max,
                      params->n_steps_field, 
                      params->dt, 
                      params->molecule.B, 
                      params->field.fwhm, 
                      params->field.e_field_squared,
                      params->molecule.delta_alpha, 
                      params->molecule.E_rot,
                      psi0,
                      params->t_start,
                      &counter, 
                      cos2,
                      weight,
                      &(params->field));

    // Solve the field-free evolution 
    field_free_propagation(params->molecule.j_max,
                           params->dt,
                           params->n_steps_field_free, 
                           &counter,
                           params->molecule.B,
                           psi0,
                           cos2,
                           weight);
}



size_t get_planar_thermal_weights(solver_params *params, double weights[params->dim])
{
    // First calculate the partition function
    size_t j_max = params->molecule.j_max;
    int j_counter = (int) -j_max;
    double partition_function = 0.0;
    for (size_t i = 0; i < params->dim; i++)
    {
        double abundance = (j_counter % 2 == 0) ? params->molecule.even_abundance : params->molecule.odd_abundance;
        partition_function += abundance * exp(-params->beta * params->molecule.E_rot[i]);
        j_counter++;
    }

    // Then find out how many states we need until the desired accuracy...
    double state_sum = 0.0;
    size_t n_states_thermal = 0;
    j_counter = 0;
    for (int i = params->molecule.j_max; i < params->dim; i++)
    {
        // ...Break if we have reached 99.9 percent
        if (state_sum / partition_function > 0.999) 
            break;

        double abundance = (double) (j_counter % 2 == 0) ? params->molecule.even_abundance : params->molecule.odd_abundance;
        double boltzmann_factor = ((i == j_max) ? 1.0 : 2.0) * abundance * exp(-params->beta * params->molecule.E_rot[i]);
        weights[i] += boltzmann_factor / partition_function;
        state_sum += boltzmann_factor;
        n_states_thermal += (i == j_max) ? 1 : 2;
        j_counter++;
    }

    return n_states_thermal;
}
