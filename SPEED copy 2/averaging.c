#include <stdio.h>
#include <math.h>
#include "averaging.h"
#include "utils.h"
#include "units.h"

#define TWO_SQRT_TWO_LOG_TWO 2.35482004503094938202

double gaussian_profile(const double r, const double I0, const double w)
{
    return I0 * exp(-2.0 * pow(r / w, 2));
}

void get_focal_average_weights(size_t n_focal, const double I0, const double w_pump, const double w_probe, double weights[n_focal], double intensities[n_focal])
{
    double r_max = 1.7 * w_probe;           // Maximum radius (Areg used 1.7? But 1.5 is more than 99% of the power)
    double dr = (double) r_max / n_focal;   // The radial step size
    double radii[n_focal];                  // The radii used form [dr; 1.5w+dr]
    double normalization = 0.0;             // The total sum of all unnormalized weights

    if (n_focal == 1)
    {
        intensities[0] = I0;
        weights[0] = 1.0;
    }
    else
    {
        for (size_t i = 0; i < n_focal; i++)
        {
            // Add 1 to avoid weight[0] = 0
            radii[i] = (double) i * dr + 0.1;
            intensities[i] = gaussian_profile(radii[i], I0, w_pump);
            weights[i] = radii[i] * gaussian_profile(radii[i], 1.0, w_probe);
            normalization += weights[i];
        }

        // Normalize the weights so they sum to unit
        for(size_t i = 0; i < n_focal; i++) 
            weights[i] /= normalization;
    }

}


double elliptical_gaussian_profile(const double I0, const double r, const double theta, const double w_x, const double w_y)
{
    return I0 * exp(-2.0 * r * r * (pow(cos(theta) / w_x, 2) + pow(sin(theta) / w_y, 2)));
}


void get_elliptical_focal_average_weights(const size_t n_r, const size_t n_theta, const double I0, const double w_x, const double w_y, 
                                            const double w_probe, double weights[n_r*n_theta], double intensities[n_r*n_theta])
{
    const double r_max = 1.7 * w_probe;
    const double dr = (double) r_max / n_r;
    const double dtheta = (double) (M_PI / 4.0) / (n_theta - 1.0);
    double normalization = 0.0;
    const double epsilon = 0.1;

    if (n_r * n_theta <= 1)
    {
        intensities[0] = I0;
        weights[0] = 1.0;
    }
    else
    {
        size_t idx = 0;
        for (size_t i = 0; i < n_r; i++)
        {
            double r = epsilon + i * dr; 
            for (int j = 0; j < n_theta; j++)
            {
                double theta = j * dtheta;
                double symmetry_factor = (theta == 0.0 || theta == M_PI / 4.0) ? 2.0 : 4.0; 
                intensities[idx] = elliptical_gaussian_profile(I0, r, theta, w_x, w_y);
                weights[idx] = symmetry_factor * r * gaussian_profile(r, 1.0, w_probe);
                normalization += weights[idx];
                idx++;
            }
        }

        for (size_t i = 0; i < n_r * n_theta; i++)
            weights[i] /= normalization;
    }
}


double gaussian_lineshape(const double B, const double B0, const double Delta_B)
{
    double sigma_B = Delta_B / TWO_SQRT_TWO_LOG_TWO;
    return exp(-0.5 * pow((B - B0) / sigma_B, 2)) / (sqrt(2.0 * M_PI) * sigma_B);
}

void get_rot_const_avg_weights(const size_t n, const double Delta_B, const double B, double weights[2*n+1], double B_values[2*n+1])
{
    // If no linewidth average! Otherwise it underestimates B
    if (n == 0)
    {
        B_values[0] = B;
        weights[0] = 1.0;
        return;
    }

    // Find the B-values
    size_t n_B = 2 * n + 1;
    double delta = 3.0 * Delta_B / n_B;
    for (size_t i = 0; i < n_B; i++)
        B_values[i] = B - 3.0 * Delta_B / 2.0 + delta * i;

    // Find the weights
    for (size_t i = 0; i < n_B; i++)
        weights[i] = gaussian_lineshape(B_values[i], B, Delta_B);
    
    // Then normalize so they sum to unity
    double weight_sum = 0.0;
    for (size_t i = 0; i < n_B; i++) weight_sum += weights[i];
    for (size_t i = 0; i < n_B; i++) weights[i] /= weight_sum;
}

double e_field_from_intensity(const double intensity)
{
    return sqrt(2.0 * intensity / (speed_of_light_SI * vaccuum_permitivity_SI)) * au_e_field_strength_per_au;
}