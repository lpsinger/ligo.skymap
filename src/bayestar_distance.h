/*
 * Copyright (C) 2015-2024  Leo Singer
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#ifndef BAYESTAR_DISTANCE_H
#define BAYESTAR_DISTANCE_H

#ifndef __cplusplus

double bayestar_distance_conditional_pdf(
    double r, double mu, double sigma, double norm);

double bayestar_distance_conditional_cdf(
    double r, double mu, double sigma, double norm);

double bayestar_distance_conditional_ppf(
    double p, double mu, double sigma, double norm);

int bayestar_distance_moments_to_parameters(
    double mean, double std, double *mu, double *sigma, double *norm);

void bayestar_distance_parameters_to_moments(
    double mu, double sigma, double *mean, double *std, double *norm);

double bayestar_volume_render(
    double x, double y, double max_distance, int axis0, int axis1,
    const double *R, long long nside, long int n, const long int *nest,
    const double *probdensity, const double *mu,
    const double *sigma, const double *norm);

double bayestar_distance_marginal_pdf(
    double r, long long npix,
    const double *prob, const double *mu,
    const double *sigma, const double *norm);

double bayestar_distance_marginal_cdf(
    double r, long long npix,
    const double *prob, const double *mu,
    const double *sigma, const double *norm);

double bayestar_distance_marginal_ppf(
    double p, long long npix,
    const double *prob, const double *mu,
    const double *sigma, const double *norm);

#endif /* __cplusplus */

#endif /* BAYESTAR_DISTANCE_H */
