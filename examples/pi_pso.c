/*
    sample_pi_pso.c -- trivial example using fpso

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
 * Sample problem.
 */

#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <limits.h>
#include "fgen.h"

/* Prototypes of characteristic functions defined for the problem. */

static double ProblemCalculateError(const FpsoPopulation *, const double *);
static void ProblemPrintIndividual(const FpsoPopulation *, const double *);

void ProblemGenerationCallback(FpsoPopulation *pop, int generation) {
	printf("Generation = %d, best error = %lf, solution ", generation, fpso_get_best_known_error(pop));
	ProblemPrintIndividual(pop, fpso_get_best_known_position(pop));
}

int main(int argc, char **argv) {
	FpsoPopulation *pop;
	pop = fpso_create(
		32,			/* Population size. */
		1,
		ProblemGenerationCallback,
		ProblemCalculateError);
	fpso_set_parameters(
		pop,
		FPSO_TOPOLOGY_GBEST,
		FPSO_BOUND_POSITION_AND_VELOCITY,
		FPSO_DEFAULT_OMEGA,
		FPSO_DEFAULT_PHI1,
		FPSO_DEFAULT_PHI2);
	fpso_set_parameter_bounds(pop, 0, 0.0, 10.0);
	fgen_random_seed_with_timer(fpso_get_rng(pop));
	fpso_run(pop, 100);
	fpso_destroy(pop);
    	exit(0);
}


/* Problem definition. */

/*
 * Example problem:
 *
 * Find PI.
 */

double ProblemCalculateError(const FpsoPopulation *pop, const double *param) {
    return fabs(param[0] - M_PI);
}

void ProblemPrintIndividual(const FpsoPopulation *pop, const double *param) {
    printf("x = %lf\n", param[0]);
}

