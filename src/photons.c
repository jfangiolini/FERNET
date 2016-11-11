/***********************************************************************************
 *                                                                                 *
 * FERNET: Fluorescence Emission Recipes and NumErical routines Toolkit            *
 * Developed by Juan F. Angiolini and Esteban Mocskos                              *
 * Facultad de Ciencias Exactas y Naturales, UBA                                   *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or (at your option) any later version.                          *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 *                                                                                 *
 ***********************************************************************************/

#include "fernet.h"

/***********************************************************************************
 * Checks if a molecule emits photons
 ***********************************************************************************/

int gaussPSF(double x, double y, double z, double w_xy, double w_z,
	     double sx, double sy, double sz, int nevents, double q, gsl_rng * r)
{
	double g, prob_abs, prob_emit;
	int phot = 0;

	g = exp(-2 * ((x - sx) * (x - sx) + (y - sy) * (y - sy)) /
		(w_xy * w_xy) - 2 * ((z - sz) * (z - sz)) / (w_z * w_z));

	for (int i = 0; i < nevents; i++) {
		prob_abs = gsl_rng_uniform(r);
		prob_emit = gsl_rng_uniform(r);
		if (g > prob_abs && prob_emit < q) {
			phot++;
		}
	}
	return phot;
}

int spimPSF(double z, double w_z, double sz, int nevents, double q, gsl_rng * r)
{
	double g, prob_abs, prob_emit;
	int phot = 0;
	g = exp(-2 * ((z - sz) * (z - sz)) / (w_z * w_z));

	for (int i = 0; i < nevents; i++) {
		prob_abs = gsl_rng_uniform(r);
		prob_emit = gsl_rng_uniform(r);
		if (g > prob_abs && prob_emit < q) {
			phot++;
		}
	}
	return phot;
}

int noiseGenerator(int photons, int noise_status, gsl_rng * r)
{
	int noise = 0;
	if (noise_status) {
		noise += gsl_ran_poisson(r, sqrt((double)photons)) + gsl_ran_gaussian_tail(r, 0, 1);
	}
	return noise;
}
