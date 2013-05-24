/*
    sample_ffit_rosenbrock.c -- ffit example for the Rosenbrock function.

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
 * Rosenbrock function example using ffit.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "fgen.h"

#define DOMAIN_X_MIN - 5.0
#define DOMAIN_X_MAX 5.0
#define DOMAIN_Y_MIN - 5.0
#define DOMAIN_Y_MAX 5.0

/* Characteristic functions defined for the Rosenbrock problem. */

static void RosenbrockGenerationCallback(Ffit *fit, int generation, const double *best_param, double best_error);
static double RosenbrockCalculateError(const Ffit *fit, const double *);

int main(int argc, char **argv) {
	Ffit *fit;
	fit = ffit_create(
		2,	/* Number of parameters. */
		RosenbrockGenerationCallback,
		RosenbrockCalculateError);
	ffit_set_parameter_range_and_mapping(fit, 0, DOMAIN_X_MIN, DOMAIN_X_MAX, FFIT_MAPPING_LINEAR);
	ffit_set_parameter_range_and_mapping(fit, 1, DOMAIN_Y_MIN, DOMAIN_Y_MAX, FFIT_MAPPING_LINEAR);
	printf("Running PSO.\n");
	ffit_run_fpso(
		fit,
		32,	/* Number of particles. */
		FPSO_TOPOLOGY_GBEST,
		FPSO_BOUND_POSITION_AND_VELOCITY,
		FPSO_DEFAULT_OMEGA,
		FPSO_DEFAULT_PHI1,
		FPSO_DEFAULT_PHI2);
	printf("Running GA.\n");
//	ffit_set_threading(fit, FFIT_THREADING_ENABLED);
	ffit_run_fgen(
		fit,
		1024,	/* Population size. */
		64,	/* Number of bits per parameter for GA, undefined for PSO. */
		FGEN_ELITIST_RANK,	/* Rank-based selection helps a lot for this problem. */
		fgen_crossover_uniform_per_bit,
		0.9,	/* Crossover rate. */
		0.015,	/* Mutation rate per bit. */
		0.005);	/* Macro-mutation rate. */
	ffit_destroy(fit);
    	exit(0);
}


/* Rosenbrock definition. */

/*
 * Example problem:
 *
 * Minimize the Rosenbrock function f(x, y) = (1 - x) ^ 2 + 100 * (y - x ^ 2) ^ 2.
 */

#define sqr(x) ((x) * (x))

static double RosenbrockCalculateError(const Ffit *fit, const double *param) {
    return sqr(1 - param[0]) + 100 * sqr(param[1] - sqr(param[0]));
}

static void RosenbrockPrintIndividual(const double *param) {
    printf("(x, y) = (%lf, %lf)\n", param[0], param[1]);
}

static void RosenbrockGenerationCallback(Ffit *fit, int generation, const double *best_param, double best_error) {
	if (generation % 50 == 0) {
		printf("Generation = %d, best error = %lf, solution ", generation, best_error);
		RosenbrockPrintIndividual(best_param);
	}
	if (generation == 500)
		ffit_signal_stop(fit);
}

