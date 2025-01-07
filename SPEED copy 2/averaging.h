#ifndef HAVE_AVERAGING_H
#define HAVE_AVERAGING_H

#include <stdio.h>

void get_focal_average_weights(size_t n_focal, const double I0, const double w_pump, const double w_probe, double weights[n_focal], double intensities[n_focal]);

void get_rot_const_avg_weights(const size_t n, const double Delta_B, const double B, double weights[2*n+1], double B_values[2*n+1]);

void get_elliptical_focal_average_weights(const size_t n_r, const size_t n_theta, const double I0, const double w_x, const double w_y, 
                                            const double w_probe, double weights[n_r*n_theta], double intensities[n_r*n_theta]);

#endif // HAVE_AVERAGING_H