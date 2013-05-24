/*
    sampe_ffit_rastrigin.c -- ffit example for a multidimensional Rastrigin function.

    fgen -- Library for optimization using a genetic algorithm or particle swarm optimization.
    Copyright 2012, Harm Hanemaaijer

    This file is part of fgen.

    fgen is free software: you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    fgen is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
    License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with fgen.  If not, see <http://www.gnu.org/licenses/>.

*/


/*
 * Rastrigin function example using ffit.
 */

#include <stdlib.h>
#include <stdio.h>
#ifndef __GNUC__
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <limits.h>
#include "fgen.h"

#define NU_DIMENSIONS 20
#define DOMAIN_MIN - 5.12
#define DOMAIN_MAX 5.12

/* Characteristic functions defined for the Rastrigin problem. */

static void RastriginGenerationCallback(Ffit *fit, int generation, const double *best_param, double best_error);
static double RastriginCalculateError(const Ffit *fit, const double *);

Ffit *setup_ffit() {
	Ffit *fit = ffit_create(
		NU_DIMENSIONS,	/* Number of parameters. */
		RastriginGenerationCallback,
		RastriginCalculateError);
	for (int i = 0; i < NU_DIMENSIONS; i++)
		ffit_set_parameter_range_and_mapping(fit, i, DOMAIN_MIN, DOMAIN_MAX, FFIT_MAPPING_LINEAR);
	return fit;
}

int main(int argc, char **argv) {
	Ffit *fit;
	fit = setup_ffit();

	printf("Running PSO.\n");
	ffit_run_fpso(
		fit,
		32,	/* Number of particles. */
		FPSO_TOPOLOGY_GBEST,
		FPSO_BOUND_POSITION_AND_VELOCITY,
		FPSO_DEFAULT_OMEGA,
		FPSO_DEFAULT_PHI1,
		FPSO_DEFAULT_PHI2);
	ffit_destroy(fit);
	fit = setup_ffit();
	ffit_set_threading(fit, FFIT_THREADING_ENABLED);
	ffit_set_generation_callback_interval(fit, 200);
	printf("Running GA (real-valued representation internally).\n");
	ffit_run_fgen_real_valued(
		fit,
		1024,
		FGEN_ELITIST_SUS,
		0.8,	/* Crossover probability. */
		0.03,	/* Mutation probability per data element. */
		0.03);	/* Macro-mutation probability per data element. */
	ffit_destroy(fit);
	fit = setup_ffit();
	printf("Running GA (bitstring representation internally).\n");
	ffit_run_fgen(
		fit,
		1024,			/* Population size. */
		64,			/* Number of bits per parameter for GA. */
		FGEN_ELITIST_SUS,	/* Rank-based selection doesn't work well with this problem. */
		fgen_crossover_uniform_per_element,
		0.9,	/* Crossover rate. */
		0.015,	/* Mutation rate per bit. */
		0.005);	/* Macro-mutation rate. */
	ffit_destroy(fit);
    	exit(0);
}


/* Problem definition. */

/*
 * Example problem:
 *
 * Minimize the n-dimensional Rastrigin function:
 *
 * A * n + SIGMA[i = 1 to n](xi ^ 2 - A * cos(2 * PI * xi))
 *
 * where A = 10.
 *
 * The two-dimensional version is given by:
 *
 * f(x, y) = 10 * 2 + x ^ 2 + y ^ 2 - 10 * cos(2 * PI * x) - 10 * cos(2 * PI * y)
 */

#define sqr(x) ((x) * (x))

static double RastriginCalculateError(const Ffit *fit, const double *param) {
	double result;
	int i;
	const double A = 10.0;
	const double n = NU_DIMENSIONS;
	result = A * n;
	for (i = 0; i < n; i++)
		result += sqr(param[i]) - A * cos(2 * M_PI * param[i]);
	return result;
}

static void RastriginPrintIndividual(const double *param) {
	int i;
	printf("(");
	for (i = 0; i < NU_DIMENSIONS - 1; i++)
		printf("%lf, ", param[i]);
	printf("%lf)\n", param[i]);
}

static void RastriginGenerationCallback(Ffit *fit, int generation, const double *best_param, double best_error) {
	if (generation % 200 == 0) {
		printf("Generation = %d, best error = %lf, solution ", generation, best_error);
		RastriginPrintIndividual(best_param);
	}
	if (generation == 1000)
		ffit_signal_stop(fit);
}

