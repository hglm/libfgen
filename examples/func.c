/*
    sample_func.c -- trivial example using fgen.

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
#include <math.h>
#include <limits.h>
#include "fgen.h"

/* Prototypes of characteristic functions defined for the problem. */

/*
 * This function returns the non-zero fitness value of a bitstring.
 */

static double ProblemCalculateFitness(const FgenPopulation *, const unsigned char *);

/*
 * Print an individual in its natural representation (e.g. not a bitstring).
 */

static void ProblemPrintIndividual(const FgenPopulation *, const FgenIndividual *);


static void ProblemGenerationCallback(FgenPopulation *pop, int generation) {
    	FgenIndividual *best = fgen_best_individual_of_population(pop);
	printf("Generation = %d, solution ", generation);
	ProblemPrintIndividual(pop, best);
}

static double FunctionUpperBound();

static double f_max;

int main(int argc, char **argv) {
	FgenPopulation *pop;
	pop = fgen_create(
		128,			/* Population size. */
		32,			/* Individual size in bits. */
		32,			/* Data element size. */
		ProblemGenerationCallback,
		ProblemCalculateFitness,
		fgen_seed_random,
		fgen_mutation_per_bit,
		fgen_crossover_one_point_per_bit);
	fgen_set_parameters(
		pop,
		FGEN_ELITIST_SUS,		/* Selection type. */
		FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
		0.600,				/* Crossover probability. */
		0.010,				/* Mutation probability per bit. */
		0);				/* Macro-mutation probability. */
	fgen_random_seed_with_timer(fgen_get_rng(pop));

	/* Calculate an upper limit of the optimum for fitness scaling. */
	f_max = FunctionUpperBound();

	fgen_run(pop, 100);
	fgen_destroy(pop);
    	exit(0);
}


/*
 * Example problem:
 *
 * Find optimum of a1 * x + a2 * x^2 .. + an * x^n.
 * in the real domain [DOMAIN_MIN, DOMAIN_MAX[.
 * A 32-bit encoding is used.
 * Resolution is (DOMAIN_MAX - DOMAIN_MIN) / 2 ^ 32.
 */

#define NU_COEFFICIENTS 2
#define DOMAIN_MIN	0.00
#define DOMAIN_MAX	10.00
#define USE_GRAY_CODING 1

static double a[NU_COEFFICIENTS] = {
	10,		/* 5x - 2x^2 */
	-2,
};

static double FunctionUpperBound() {
	double f;
	double dom_max;
	double x_power;
	int i;
	dom_max = fabs(DOMAIN_MIN);
	if (fabs(DOMAIN_MAX) > dom_max)
		dom_max = fabs(DOMAIN_MAX);
	f = 0;
	/* x = dommax */
	x_power = dom_max;
	for (i = 0; i < NU_COEFFICIENTS; i++) {
		f += x_power * fabs(a[i]);
		x_power *= dom_max;	/* xpower = x ^ (i + 1) */
	}
/*        printf("Upper Bound = %lf\n", f); */
	return f;
}

static double FunctionValue( double x ) {
	double f;
	double x_power;
	int i;
	/* Calculate the function value. */
	f = 0;
	x_power = x;
	for (i = 0; i < NU_COEFFICIENTS; i++) {
		f += x_power * a[i];
		x_power *= x;		/* xpower = x ^ (i + 1) */
	}
	return f;
}


static double ProblemCalculateFitness(const FgenPopulation *pop, const unsigned char *bitstring) {
    double x, f;
    if (USE_GRAY_CODING) {
        unsigned char decoded[8];
        fgen_decode_from_gray(pop, bitstring, decoded);
        x = fgen_bitstring_uint32_to_double(decoded, DOMAIN_MIN, DOMAIN_MAX);
    }
    else
        x = fgen_bitstring_uint32_to_double(bitstring, DOMAIN_MIN, DOMAIN_MAX);

    f = FunctionValue(x);
    /* Scale such that negative numbers have lower fitness. */
    return (f + f_max) / (f_max * 2);
}

static void ProblemPrintIndividual(const FgenPopulation *pop, const FgenIndividual *ind) {
    double x;
    if (USE_GRAY_CODING) {
        unsigned char decoded[4];
        fgen_decode_from_gray(pop, ind->bitstring, decoded);
        x = fgen_bitstring_uint32_to_double(decoded, DOMAIN_MIN, DOMAIN_MAX);
    }
    else
        x = fgen_bitstring_uint32_to_double(ind->bitstring, DOMAIN_MIN, DOMAIN_MAX);
    printf("x = %lf, f(x) = %lf\n", x, FunctionValue(x));
}

