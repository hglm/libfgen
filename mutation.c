
/*
    mutation.c -- Genetic algorithm mutation implementation.

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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include "fgen.h"		/* General includes. */
#include "parameters.h"
#include "population.h"		/* Module includes. */
#include "bitstring.h"
#include "random.h"
#include "win32_compat.h"

#ifndef __GNUC__

// Log-gamma code was taken from the following source:
// Visit http://www.johndcook.com/stand_alone_code.html for the source of this code and more like it.

double LogGamma(double x);

static double Gamma
(
    double x    // We require x > 0
)
{

    // Split the function domain into three intervals:
    // (0, 0.001), [0.001, 12), and (12, infinity)

    ///////////////////////////////////////////////////////////////////////////
    // First interval: (0, 0.001)
	//
	// For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
	// So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
	// The relative error over this interval is less than 6e-7.

	const double gamma = 0.577215664901532860606512090; // Euler's gamma constant

    if (x < 0.001)
        return 1.0/(x*(1.0 + gamma*x));

    ///////////////////////////////////////////////////////////////////////////
    // Second interval: [0.001, 12)
    
	if (x < 12.0)
    {
        // The algorithm directly approximates gamma over (1,2) and uses
        // reduction identities to reduce other arguments to this interval.
		
		double y = x;
        int n = 0;
        bool arg_was_less_than_one = (y < 1.0);

        // Add or subtract integers as necessary to bring y into (1,2)
        // Will correct for this below
        if (arg_was_less_than_one)
        {
            y += 1.0;
        }
        else
        {
            n = static_cast<int> (floor(y)) - 1;  // will use n later
            y -= n;
        }

        // numerator coefficients for approximation over the interval (1,2)
        static const double p[] =
        {
            -1.71618513886549492533811E+0,
             2.47656508055759199108314E+1,
            -3.79804256470945635097577E+2,
             6.29331155312818442661052E+2,
             8.66966202790413211295064E+2,
            -3.14512729688483675254357E+4,
            -3.61444134186911729807069E+4,
             6.64561438202405440627855E+4
        };

        // denominator coefficients for approximation over the interval (1,2)
        static const double q[] =
        {
            -3.08402300119738975254353E+1,
             3.15350626979604161529144E+2,
            -1.01515636749021914166146E+3,
            -3.10777167157231109440444E+3,
             2.25381184209801510330112E+4,
             4.75584627752788110767815E+3,
            -1.34659959864969306392456E+5,
            -1.15132259675553483497211E+5
        };

        double num = 0.0;
        double den = 1.0;
        int i;

        double z = y - 1;
        for (i = 0; i < 8; i++)
        {
            num = (num + p[i])*z;
            den = den*z + q[i];
        }
        double result = num/den + 1.0;

        // Apply correction if argument was not initially in (1,2)
        if (arg_was_less_than_one)
        {
            // Use identity gamma(z) = gamma(z+1)/z
            // The variable "result" now holds gamma of the original y + 1
            // Thus we use y-1 to get back the orginal y.
            result /= (y-1.0);
        }
        else
        {
            // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
            for (i = 0; i < n; i++)
                result *= y++;
        }

		return result;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Third interval: [12, infinity)

    if (x > 171.624)
    {
		// Correct answer too large to display. Force +infinity.
		double temp = DBL_MAX;
//		return temp*2.0;
		return temp;
    }

    return exp(LogGamma(x));
}

static double LogGamma
(
    double x    // x must be positive
)
{

    if (x < 12.0)
    {
        return log(fabs(Gamma(x)));
    }

	// Abramowitz and Stegun 6.1.41
    // Asymptotic series should be good to at least 11 or 12 figures
    // For error analysis, see Whittiker and Watson
    // A Course in Modern Analysis (1927), page 252

    static const double c[8] =
    {
		 1.0/12.0,
		-1.0/360.0,
		1.0/1260.0,
		-1.0/1680.0,
		1.0/1188.0,
		-691.0/360360.0,
		1.0/156.0,
		-3617.0/122400.0
    };
    double z = 1.0/(x*x);
    double sum = c[7];
    for (int i=6; i >= 0; i--)
    {
        sum *= z;
        sum += c[i];
    }
    double series = sum/x;

    static const double halfLogTwoPi = 0.91893853320467274178032973640562;
    double logGamma = (x - 0.5)*log(x) - x + halfLogTwoPi + series;    
	return logGamma;
}

#define lgamma(x) LogGamma(x)

#endif

/* Avoid function calls to fgen_mutate_bit. */
#define mutate_bit(bitstring, bitnumber) bitstring[bitnumber / 8] ^= 1 << (bitnumber & 7);

static double BinomialCoefficient(int n, int k) {
	return exp(lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1));
}

/*
 * This function is called at the beginning of fgen_run and variants. It sets up an array of probabilities of
 * how many bits in an individual will be mutated.
 */

void SetupFastMutationCumulativeChanceArray(FgenPopulation *pop) {
	if (pop->fast_mutation_cumulative_chance != NULL) {
		// There is already a table present.
		if (pop->mutation_probability_float == pop->fast_mutation_probability)
			// We have already calculated the table for the right probability.
			return;
		// Otherwise, overwrite the existing table.
	}
	else
		// There was no table yet.
		pop->fast_mutation_cumulative_chance = (float *)malloc(pop->individual_size_in_bits * sizeof(float));
	/* Calculate the chance that one or more bits are mutated in an individual. Equal to one minus the */
	/* minus the chance that no bits are mutated. */
	pop->fast_mutation_cumulative_chance[0] =
		1.0 - pow(1.0 - pop->mutation_probability_float, pop->individual_size_in_bits);
//	printf("Probability %d bits or more mutated: %lf\n", 1, pop->fast_mutation_cumulative_chance[0]);
	for (int i = 1; i < pop->individual_size_in_bits; i++) {
		/* Calculate the chance that exactly i bits are mutated in an individual. */
		float f = BinomialCoefficient(pop->individual_size_in_bits, i) *
			pow(pop->mutation_probability_float, i) *
			pow(1.0 - pop->mutation_probability_float, pop->individual_size_in_bits - i);
		if (isnan(f))
			pop->fast_mutation_cumulative_chance[i] = 0;
		else
			/* Assign the chance that more than i bits are mutated. */
			pop->fast_mutation_cumulative_chance[i] = pop->fast_mutation_cumulative_chance[i - 1] - f;
//		printf("Probability %d bits or more mutated: %lf\n", i + 1, pop->fast_mutation_cumulative_chance[i]);
		if (pop->fast_mutation_cumulative_chance[i] <= 0.000001) {
			if (i < pop->individual_size_in_bits - 1)
				pop->fast_mutation_cumulative_chance[i + 1]  = 0;
//			printf("Probability %d bits or more mutated: %lf\n", i + 2,
//				pop->fast_mutation_cumulative_chance[i + 1]);
			break;
		}
	}
	pop->fast_mutation_probability = pop->mutation_probability_float;
}

static int NumberOfBitsToMutate(FgenPopulation *pop) {
	float f = fgen_random_f(pop->rng, 1.0);
	int n = 0;
	while (f < pop->fast_mutation_cumulative_chance[n])
		n++;
	return n;
}

/*
 * This function performs the Mutation step of the genetic algorithm.
 */

void DoMutation(FgenPopulation *pop) {
	int i;
	FgenIndividual **ind1 = pop->ind;
	FgenIndividual **ind2 = AllocatePopulation(pop);
	if (pop->fgen_mutation_func == fgen_mutation_per_bit_fast ||
	pop->fgen_mutation_func == fgen_mutation_per_bit_plus_macro_mutation_fast) {
		for (i = 0; i < pop->size; i++) {
			if ((pop->selection_type & FGEN_ELITIST_ELEMENT) && ind1[i]->is_elite) {
				ind2[i] = ind1[i];
				/* The following two lines of code have no net effect so can be omitted. */
				/* ind1[i]->refcount++; */
				/* FreeIndividual(ind1[i]); */
				continue;
			}
			int n = NumberOfBitsToMutate(pop);
			if (n == 0 && pop->fgen_mutation_func == fgen_mutation_per_bit_fast)
				ind2[i] = ind1[i];
			else {
				ind2[i] = NewIndividual(pop);
				CopyIndividualBitstring(pop, ind1[i], ind2[i]);
				pop->fast_mutation_nu_bits_to_mutate = n;
				pop->fgen_mutation_func(pop, ind1[i]->bitstring, ind2[i]->bitstring);
				FreeIndividual(ind1[i]);
			}
		}
	}
	else {
		for (i = 0; i < pop->size; i++) {
			if ((pop->selection_type & FGEN_ELITIST_ELEMENT) && ind1[i]->is_elite) {
				ind2[i] = ind1[i];
				/* The following two lines of code have no net effect so can be omitted. */
				/* ind1[i]->refcount++; */
				/* FreeIndividual(ind1[i]); */
				continue;
			}
			ind2[i] = NewIndividual(pop);
			CopyIndividualBitstring(pop, ind1[i], ind2[i]);
			pop->fgen_mutation_func(pop, ind1[i]->bitstring, ind2[i]->bitstring);
			FreeIndividual(ind1[i]);
		}
	}
	free(ind1);
	pop->ind = ind2;
}

/**
 * Bit-wise mutation operator. Each bit is toggled with the preset mutation probability per bit.
 */

void fgen_mutation_per_bit(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	int i;
	for (i = 0; i < pop->individual_size_in_bits; i++) {
		int r = RandomBits(pop->rng, 16);
		if (r < pop->mutation_probability) {
			/* Mutate the bit. */
			mutate_bit(child, i);
		}
	}
}

/**
 * Fast bit-wise mutation operator. A predefined random number of bits is mutated. There is a chance that the same bit
 * is toggled more than once, but because the mutation rate per bit is low it is not that important.
 * pop->fast_mutation_nu_bits_to_mutate is guaranteed to be >= 1.
 */

void fgen_mutation_per_bit_fast(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	int n = pop->fast_mutation_nu_bits_to_mutate;
	if (pop->individual_size_shift >= 0) {
		// Optimize for individual size in bits that is a power of two.
		// Individual size in bits is (1 << pop->individual_size_shift).
		unsigned int r = RandomIntPowerOfTwoWithShift(pop->rng,
			pop->individual_size_shift);
		mutate_bit(child, r);
                for (int i = 1; i < n; i++) {
			r = RandomIntPowerOfTwoRepeat(pop->rng);
			mutate_bit(child, r);
		}
		return;
        }
        RandomIntGeneralPrepareForRepeat(pop->rng, pop->individual_size_in_bits);
	for (int i = 0; i < n; i++) {
		int r = RandomIntGeneralRepeat(pop->rng);
		mutate_bit(child, r);
	}
}

/**
 * Mutation operator that performs both a bit-wise mutation and a macro mutation where a random data element
 * is entirely replaced with a random value. The preset macro mutation probability is applied per data element.
 */

void fgen_mutation_per_bit_plus_macro_mutation(FgenPopulation *pop, const unsigned char *parent,
unsigned char *child) {
	fgen_mutation_per_bit(pop, parent, child);
	if (pop->macro_mutation_probability > 0) {
		int i;
		/* For each data element. */
		int n = pop->individual_size_in_bits / pop->data_element_size;
		for (i = 0; i < n; i++) {
			int r = RandomBits(pop->rng, 16);
			if (r < pop->macro_mutation_probability)
				/* Mutate the entire data element with a random value */
				CreateRandomBitstring(pop->rng, child + i * pop->data_element_size / 8,
					pop->data_element_size / 8);
		}
	}
}

/**
 * Mutation operator that performs both a bit-wise mutation and a macro mutation where a random data element
 * is entirely replaced with a random value. The preset macro mutation probability is applied per data element.
 * Fast version.
 */

void fgen_mutation_per_bit_plus_macro_mutation_fast(FgenPopulation *pop, const unsigned char *parent,
unsigned char *child) {
	fgen_mutation_per_bit_fast(pop, parent, child);
	if (pop->macro_mutation_probability > 0) {
		if (pop->data_element_size == 32) {
			int n = pop->individual_size_in_bits / 32;
			for (int i = 0; i < n; i++) {
				int r = RandomBits(pop->rng, 16);
				if (r < pop->macro_mutation_probability)
					*(unsigned int *)&child[i * 4] = fgen_random_32(pop->rng);
			}
			return;
		}
		if (pop->data_element_size == 64) {
			int n = pop->individual_size_in_bits / 32;
			for (int i = 0; i < n; i += 2) {
				int r = RandomBits(pop->rng, 16);
				if (r < pop->macro_mutation_probability) {
					*(unsigned int *)&child[i * 4] = fgen_random_32(pop->rng);
					*(unsigned int *)&child[i * 4 + 4] = fgen_random_32(pop->rng);
				}
			}
			return;
		}
		/* For each data element. */
		int n = pop->individual_size_in_bits / pop->data_element_size;
		for (int i = 0; i < n; i++) {
			int r = RandomBits(pop->rng, 16);
			if (r < pop->macro_mutation_probability)
				/* Mutate the entire data element with a random value */
				CreateRandomBitstring(pop->rng, child + i * pop->data_element_size / 8,
					pop->data_element_size / 8);
		}
	}
}

/**
 * Simple mutation operator for permutations.
 */

void fgen_mutation_permutation_swap(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	int index1, index2;
	int r = RandomBits(pop->rng, 16);
	if (r >= pop->mutation_probability)
		return;
        RandomIntGeneralPrepareForRepeat(pop->rng, pop->permutation_size);
	do {
		index1 = RandomIntGeneralRepeat(pop->rng);
		index2 = RandomIntGeneralRepeat(pop->rng);
	} while (index1 == index2);
	*(int *)&child[index1 * 4] = *(int *)&parent[index2 * 4];
	*(int *)&child[index2 * 4] = *(int *)&parent[index1 * 4];
}

/**
 * Insert mutation operator for permutations. Takes a random element from the permutation and reinserts it at a random
 * location in the permutation.
 */

void fgen_mutation_permutation_insert(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	int r = RandomBits(pop->rng, 16);
	if (r >= pop->mutation_probability)
		return;
        RandomIntGeneralPrepareForRepeat(pop->rng, pop->permutation_size);
	int index = RandomIntGeneralRepeat(pop->rng);
	int index_dest = RandomIntGeneralRepeat(pop->rng);
	if (index_dest < index) {
		/* The destination of the element is to the left of its previous position. */
		/* The first, unaltered part of the permutation was already copied into the child. */
		/* The second part is the inserted element. */
		*(int *)&child[index_dest * 4] = *(int *)&parent[index * 4];
		/* The third part is the part of the parent up to the element location. */
		memcpy(child + (index_dest + 1) * 4, parent + index_dest * 4, (index - index_dest) * 4);
		/* The fourth, unaltered part of the permutation was already copied into the child. */
	}
	if (index_dest > index) {
		/* The destination of the element is to the right of it previous position. */
		/* The first, unaltered part of the permutation was already copied into the child. */
		/* The second part is the part of the parent to the right of the element, up to destination location. */
		memcpy(child + index * 4, parent + (index + 1) * 4, (index_dest - index) * 4);
		/* The third part is the inserted element. */
		*(int *)&child[index_dest * 4] = *(int *)&parent[index * 4];
		/* The fourth, unaltered part of the permutation was already copied into the child. */
	}
	/* If the destination is the same as the source location, nothing happens (the child already contains a copy). */
}

/**
 * Invert mutation operator for permutations that inverts a random subroute within the permutation.
 */

void fgen_mutation_permutation_invert(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	int r = RandomBits(pop->rng, 16);
	if (r >= pop->mutation_probability)
		return;
	const int *p = (int *)parent;
	int *c = (int *)child;
	int index = RandomInt(pop->rng, pop->permutation_size);
	int max_size = pop->permutation_size - index;
	int size;
	if (max_size > 1)
		size = RandomInt(pop->rng, max_size) + 1;
	else
		size = 1;
	/* Invert the subroute. */
	for (int i = index; i < index + size; i++) {
		c[i] = p[index + size - 1 - (i - index)];
	}
}

/*
 * Mutation, threaded version.
 * It is slower, and may not be safe. Not used at the moment.
 */

typedef struct {
	FgenPopulation *pop;
	FgenIndividual **ind2;
	int start_index;
	int end_index;
} ThreadData;

void *DoMutationPartThread(void *threadarg) {
	ThreadData *thread_data = (ThreadData *)threadarg;
	FgenPopulation *pop = thread_data->pop;
	FgenIndividual **ind2 = thread_data->ind2;
	int start_index = thread_data->start_index;
	int end_index = thread_data->end_index;
	FgenIndividual **ind1 = pop->ind;
	int i;
	for (i = start_index; i < end_index; i++) {
		if ((pop->selection_type & FGEN_ELITIST_ELEMENT) && ind1[i]->is_elite) {
			ind2[i] = ind1[i];
			ind1[i]->refcount++;	/* Is this safe? */
			continue;
		}
		ind2[i] = NewIndividual(pop);
		CopyIndividualBitstring(pop, ind1[i], ind2[i]);
		pop->fgen_mutation_func(pop, ind1[i]->bitstring, ind2[i]->bitstring);
	}
	return NULL;
}

void DoMutationThreaded(FgenPopulation *pop) {
	ThreadData *thread_data;
	pthread_t *thread;
	int max_threads;
	FgenIndividual **ind2;
	int i;
	max_threads = 2;
	thread = (pthread_t *)malloc(sizeof(pthread_t) * max_threads);
	thread_data = (ThreadData *)malloc(sizeof(ThreadData) * max_threads);
	ind2 = AllocatePopulation(pop);
	/* Do the first half. */
	thread_data[0].pop = pop;
	thread_data[0].ind2 = ind2;
	thread_data[0].start_index = 0;
	thread_data[0].end_index = pop->size / 2;
	pthread_create(&thread[0], NULL, DoMutationPartThread, &thread_data[0]);
	/* Do the second half. */
	thread_data[1].pop = pop;
	thread_data[1].ind2 = ind2;
	thread_data[1].start_index = pop->size / 2;
	thread_data[1].end_index = pop->size;
	DoMutationPartThread(&thread_data[1]);
	pthread_join(thread[0], NULL);
	for (i = 0; i < pop->size; i++)
		FreeIndividual(pop->ind[i]);
	free(thread_data);
	free(thread);
	free(pop->ind);
	pop->ind = ind2;
}

