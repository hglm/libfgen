/*
    crossover.c -- Crossover stage of the genetic algorithm.

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
#include <malloc.h>	// Added for Visual C++
#include "string.h"
#include "fgen.h"
#include "parameters.h"
#include "population.h"
#include "random.h"
#include "bitstring.h"
#include "error.h"


/*
 * This function performs the Crossover step of the genetic algorithm.
 * The populationSize must be a multiple of 2.
 */

void DoCrossover(FgenPopulation *pop) {
	int i;
	if (pop->crossover_probability == 0 || pop->fgen_crossover_func == fgen_crossover_noop)
		return;
	FgenIndividual **ind1 = pop->ind;
	FgenIndividual **ind2 = AllocatePopulation(pop);
	for (i = 0; i < pop->size - 1; i += 2) {
		if ((pop->selection_type & FGEN_ELITIST_ELEMENT) &&
		(ind1[i]->is_elite || ind1[i + 1]->is_elite)) {
			ind2[i] = ind1[i];
			ind2[i + 1] = ind1[i + 1];
			/* The following four lines of code have no net effect and can be omitted. */
			/* ind1[i]->refcount++; */
			/* ind1[i + 1]->refcount++; */
			/* FreeIndividual(ind1[i]); */
			/* FreeIndividual(ind1[i + 1]); */
			continue;
		}
		int r;
		r = RandomBits(pop->rng, 8);
		if (r < pop->crossover_probability) {
			/* Perform crossover operation. */
			ind2[i] = NewIndividual(pop);
			ind2[i + 1] = NewIndividual(pop);
			pop->fgen_crossover_func(pop, ind1[i]->bitstring, ind1[i + 1]->bitstring, ind2[i]->bitstring,
				ind2[i + 1]->bitstring);
			FreeIndividual(ind1[i]);
			FreeIndividual(ind1[i + 1]);
		}
		else {
			/* Assign the unmodified individuals. */
			ind2[i] = ind1[i];
			ind2[i + 1] = ind1[i + 1];
			/* The following four lines of code have no net effect and can be omitted. */
			/* ind1[i]->refcount++; */
			/* ind1[i + 1]->refcount++; */
			/* FreeIndividual(ind1[i]); */
			/* FreeIndividual(ind1[i + 1]); */
		}
	}
	free(ind1);
	pop->ind = ind2;
}


/**
 * Per element one-point crossover operator that respects data element boundaries.
 */

void fgen_crossover_one_point_per_element(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2) {
	int bitnumber;
	bitnumber = RandomInt(pop->rng, pop->nu_data_elements) *
		pop->data_element_size;
	if (pop->data_element_size_shift >= 3) {
                // data_element_size is power of two >= 8.
		memcpy(child1, parent1, bitnumber / 8);
		memcpy(child1 + bitnumber / 8, parent2 + bitnumber / 8, (pop->individual_size_in_bits - bitnumber) / 8);
		memcpy(child2, parent2, bitnumber / 8);
		memcpy(child2 + bitnumber / 8, parent1 + bitnumber / 8, (pop->individual_size_in_bits - bitnumber) / 8);
		return;
	}
	/* It is quickest to copy the whole string first. */
	memcpy(child1, parent1, INDIVIDUAL_SIZE_IN_BYTES(pop));
	memcpy(child2, parent2, INDIVIDUAL_SIZE_IN_BYTES(pop));
	fgen_copy_partial_bitstring(parent1, child2, bitnumber,
		pop->individual_size_in_bits - bitnumber);
	fgen_copy_partial_bitstring(parent2, child1, bitnumber,
		pop->individual_size_in_bits - bitnumber);
}

/**
 * Bit-wise one-point crossover operator.
 */

void fgen_crossover_one_point_per_bit(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2) {
	int bitnumber;
	bitnumber = RandomInt(pop->rng, pop->individual_size_in_bits);
	/* It is quickest to copy the whole string first. */
	memcpy(child1, parent1, INDIVIDUAL_SIZE_IN_BYTES(pop));
	memcpy(child2, parent2, INDIVIDUAL_SIZE_IN_BYTES(pop));
	fgen_copy_partial_bitstring(parent1, child2, bitnumber,
		pop->individual_size_in_bits - bitnumber);
	fgen_copy_partial_bitstring(parent2, child1, bitnumber,
		pop->individual_size_in_bits - bitnumber);
}

/** 
 * Per-element two-point crossover operator that respects data element boundaries.
 */

void fgen_crossover_two_point_per_element(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2) {
	/* Exchange the segments that fall between two positions. */
	int bitnumber1, bitnumber2;
        RandomIntGeneralPrepareForRepeat(pop->rng, pop->nu_data_elements);
	bitnumber1 = RandomIntGeneralRepeat(pop->rng) * pop->data_element_size;
	bitnumber2 = RandomIntGeneralRepeat(pop->rng) * pop->data_element_size;
	if (bitnumber1 > bitnumber2) {
		int temp;
		temp = bitnumber1;
		bitnumber1 = bitnumber2;
		bitnumber2 = temp;
	}
	if (pop->data_element_size_shift >= 3) {
                // data_element_size is power of two >= 8.
		memcpy(child1, parent1, bitnumber1 / 8);
		memcpy(child1 + bitnumber1 / 8, parent2 + bitnumber1 / 8, (bitnumber2 - bitnumber1) / 8);
		memcpy(child1 + bitnumber2 / 8, parent1 + bitnumber2 / 8, (pop->individual_size_in_bits - bitnumber2) / 8);
		memcpy(child2, parent2, bitnumber1 / 8);
		memcpy(child2 + bitnumber1 / 8, parent1 + bitnumber1 / 8, (bitnumber2 - bitnumber1) / 8);
		memcpy(child2 + bitnumber2 / 8, parent2 + bitnumber2 / 8, (pop->individual_size_in_bits - bitnumber2) / 8);
		return;
	}
	/* It is quickest to copy the whole string first. */
	memcpy(child1, parent1, INDIVIDUAL_SIZE_IN_BYTES(pop));
	memcpy(child2, parent2, INDIVIDUAL_SIZE_IN_BYTES(pop));
	fgen_copy_partial_bitstring(parent1, child2, bitnumber1,
		bitnumber2 - bitnumber1);
	fgen_copy_partial_bitstring(parent2, child1, bitnumber1,
		bitnumber2 - bitnumber1);
}

/**
 * Bit-wise two-point crossover operator.
 */

void fgen_crossover_two_point_per_bit(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2) {
	/* Exchange the segments that fall between two positions. */
	int bitnumber1, bitnumber2;
	RandomIntGeneralPrepareForRepeat(pop->rng, pop->individual_size_in_bits);
	bitnumber1 = RandomIntGeneralRepeat(pop->rng);
	bitnumber2 = RandomIntGeneralRepeat(pop->rng);
	if (bitnumber1 > bitnumber2) {
		int temp;
		temp = bitnumber1;
		bitnumber1 = bitnumber2;
		bitnumber2 = temp;
	}
	/* It is quickest to copy the whole string first. */
	memcpy(child1, parent1, INDIVIDUAL_SIZE_IN_BYTES(pop));
	memcpy(child2, parent2, INDIVIDUAL_SIZE_IN_BYTES(pop));
	fgen_copy_partial_bitstring(parent1, child2, bitnumber1,
		bitnumber2 - bitnumber1);
	fgen_copy_partial_bitstring(parent2, child1, bitnumber1,
		bitnumber2 - bitnumber1);
}

/** 
 * Bit-wise uniform crossover operator. Each bit is randomly selected from one of the two parents.
 */

#define uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, mask, r) \
	destword1 |= (srcword1 & mask) & (r & mask); \
	destword2 |= (srcword2 & mask) & (r & mask); \
	destword1 |= (srcword2 & mask) & ((r & mask) ^ mask); \
	destword2 |= (srcword1 & mask) & ((r & mask) ^ mask);

void fgen_crossover_uniform_per_bit(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2) {
	/* Randomly select the parent for each bit. */
	/* pop->individual_size_in_bits is assumed to be multiple of 8. */
	unsigned char *srcstr1, *srcstr2;
	unsigned char *deststr1, *deststr2;
	int i;
	srcstr1 = (unsigned char *)parent1;
	srcstr2 = (unsigned char *)parent2;
	deststr1 = child1;
	deststr2 = child2;
        int n = INDIVIDUAL_SIZE_IN_BYTES(pop);
        if ((n & 3) != 0)
		goto not_multiple_of_four_bytes;
        /* Multiple of four bytes (32-bit words) allows further optimization. */
        n >>= 2;
	for (i = 0; i < n; i++) {
		unsigned int srcword1, srcword2;
		unsigned int destword1, destword2;
		srcword1 = *((unsigned int *)srcstr1 + i);
		srcword2 = *((unsigned int *)srcstr2 + i);
		destword1 = 0x00;
		destword2 = 0x00;
		unsigned int r = fgen_random_32(pop->rng);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x1, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x2, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x4, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x8, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x10, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x20, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x40, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x80, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x100, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x200, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x400, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x800, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x1000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x2000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x4000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x8000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x10000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x20000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x40000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x80000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x100000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x200000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x400000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x800000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x1000000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x2000000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x4000000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x8000000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x10000000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x20000000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x40000000, r);
		uniform_crossover_merge_bit(srcword1, srcword2, destword1, destword2, 0x80000000, r);
		*((unsigned int *)deststr1 + i) = destword1;
		*((unsigned int *)deststr2 + i) = destword2;
	}
	return;
not_multiple_of_four_bytes :
	for (i = 0; i < INDIVIDUAL_SIZE_IN_BYTES(pop); i++) {
		unsigned int srcbyte1, srcbyte2;
		unsigned int destbyte1, destbyte2;
		srcbyte1 = srcstr1[i];
		srcbyte2 = srcstr2[i];
		destbyte1 = 0x00;
		destbyte2 = 0x00;
		int r = RandomBits(pop->rng, 8);
		uniform_crossover_merge_bit(srcbyte1, srcbyte2, destbyte1, destbyte2, 0x1, r);
		uniform_crossover_merge_bit(srcbyte1, srcbyte2, destbyte1, destbyte2, 0x2, r);
		uniform_crossover_merge_bit(srcbyte1, srcbyte2, destbyte1, destbyte2, 0x4, r);
		uniform_crossover_merge_bit(srcbyte1, srcbyte2, destbyte1, destbyte2, 0x8, r);
		uniform_crossover_merge_bit(srcbyte1, srcbyte2, destbyte1, destbyte2, 0x10, r);
		uniform_crossover_merge_bit(srcbyte1, srcbyte2, destbyte1, destbyte2, 0x20, r);
		uniform_crossover_merge_bit(srcbyte1, srcbyte2, destbyte1, destbyte2, 0x40, r);
		uniform_crossover_merge_bit(srcbyte1, srcbyte2, destbyte1, destbyte2, 0x80, r);
		deststr1[i] = destbyte1;
		deststr2[i] = destbyte2;
	}
}

/**
 * Crossover operator that respects data element boundaries. Each data element is selected randomly from one of
 * the two parents. This is also called discrete recombination.
 */

void fgen_crossover_uniform_per_element(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2) {
	/* Randomly select the parent for each data element. */
	/* pop->individual_size_in_bits is assumed to be multiple of pop->data_element_size. */
	if (pop->data_element_size == 32) {
		unsigned int *p1 = (unsigned int *)parent1;
		unsigned int *p2 = (unsigned int *)parent2;
		unsigned int *c1 = (unsigned int *)child1;
		unsigned int *c2 = (unsigned int *)child2;
		for (int i = 0; i < pop->individual_size_in_bits / 32; i++) {
			int r = RandomBits(pop->rng, 1);
			if (r == 0) {
				/* Swap. */
				c1[i] = p2[i];
				c2[i] = p1[i];
			}
			else {
				c1[i] = p1[i];
				c2[i] = p2[i];
			}
		}
		return;
	}
	if (pop->data_element_size == 64) {
		unsigned int *p1 = (unsigned int *)parent1;
		unsigned int *p2 = (unsigned int *)parent2;
		unsigned int *c1 = (unsigned int *)child1;
		unsigned int *c2 = (unsigned int *)child2;
		for (int i = 0; i < pop->individual_size_in_bits / 32; i += 2) {
			int r = RandomBits(pop->rng, 1);
			if (r == 0) {
				/* Swap. */
				c1[i] = p2[i];
				c1[i + 1] = p2[i + 1];
				c2[i] = p1[i];
				c2[i + 1] = p1[i + 1];
			}
			else {
				c1[i] = p1[i];
				c1[i + 1] = p1[i + 1];
				c2[i] = p2[i];
				c2[i + 1] = p2[i + 1];
			}
		}
		return;
	}
	for (int i = 0; i < pop->individual_size_in_bits; i += pop->data_element_size) {
		int r = RandomBits(pop->rng, 1);
		if (r == 0) {
			/* Swap this element. */
			fgen_copy_partial_bitstring(parent2, child1, i, pop->data_element_size);
			fgen_copy_partial_bitstring(parent1, child2, i, pop->data_element_size);
		}
		else {
			/* Copy elements unchanged. */
			fgen_copy_partial_bitstring(parent1, child1, i, pop->data_element_size);
			fgen_copy_partial_bitstring(parent2, child2, i, pop->data_element_size);
		}
	}
}

/**
 * Order based (OX1) crossover operator for permutations.
 */

void fgen_crossover_permutation_order_based(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2) {
	const int *p1 = (int *)parent1;
	const int *p2 = (int *)parent2;
	int *c1 = (int *)child1;
	int *c2 = (int *)child2;
	// Pick a random subroute of random size.
	int index = RandomInt(pop->rng, pop->permutation_size);
	int max_size = pop->permutation_size - index;
        int size;
	if (max_size == pop->permutation_size)
		size = RandomInt(pop->rng, pop->permutation_size - 1) + 1;
	else
        if (max_size > 1)
		size = RandomInt(pop->rng, max_size) + 1;
	else
		size = 1;
	// Copy the subroute from parent 1 to the same location in child 1.
	memcpy(c1 + index, p1 + index, size * 4);
	// Mark out elements in parent 2 that exist in the subroute.
	int *skip = (int *)alloca(pop->permutation_size * sizeof(int));
	for (int i = 0; i < pop->permutation_size; i++)
		skip[i] = 0;
	for (int i = index; i < index + size; i++)
		skip[p1[i]] = 1;
	// Starting on the right side of the subroute, take elements from parent 2 and insert them at the right side
	// of the subroute in child 1, skipping elements in parent 2 that were marked out.
	int i = index + size;
	int j = index + size;
	if (j == pop->permutation_size)
		j = 0;
	for (;;) {
		if (i == pop->permutation_size)
			i = 0;
		if (!skip[p2[i]]) {
			c1[j] = p2[i];
			j++;
			if (j == pop->permutation_size)
				j = 0;
			if (j == index)
				break;
		}
		i++;
	}
	// Now do the same for the second child.
	// Copy the subroute from parent 2 to the same location in child 2.
	memcpy(c2 + index, p2 + index, size * 4);
	// Mark out elements in parent 1 that exist in the subroute.
	for (int i = 0; i < pop->permutation_size; i++)
		skip[i] = 0;
	for (int i = index; i < index + size; i++)
		skip[p2[i]] = 1;
	// Starting on the right side of the subroute, take elements from parent 1 and insert them at the right side
	// of the subroute in child 2, skipping elements in parent 1 that were marked out.
	i = index + size;
	j = index + size;
	if (j == pop->permutation_size)
		j = 0;
	for (;;) {
		if (i == pop->permutation_size)
			i = 0;
		if (!skip[p1[i]]) {
			c2[j] = p1[i];
			j++;
			if (j == pop->permutation_size)
				j = 0;
			if (j == index)
				break;
		}
		i++;
	}
}

static inline int TestNumber(const int *perm, int size, int value) {
	int i;	
	for (i = 0; i < size; i++) {
		if (perm[i] == value)
			return 1;
	}
	return 0;
}

/** 
 * Position-based crossover operator for permutations. Effective for the travelling salesman problem.
 */

void fgen_crossover_permutation_position_based(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2) {
	int i;
	const int *p1 = (int *)parent1;
	const int *p2 = (int *)parent2;
	int *c1 = (int *)child1;
	int *c2 = (int *)child2;
	int pos_c1;
	int pos_c2;
	int size = pop->permutation_size;
	for (i = 0; i < size; i++) {
		c1[i] = - 1;
		c2[i] = - 1;
	}
	/* Pick a set of random cities. */
	int pos = RandomInt(pop->rng, size - 2);
	/* Keep adding random cities until we can add no more. Copy the selected cities into the offspring. */
	while (pos < size) {
		c1[pos] = p1[pos];
		c2[pos] = p2[pos];
		pos += RandomInt(pop->rng, size - pos) + 1;
	}
	/* Fill in the blanks. First create two position markers so that we know where we are in */
	/* c1 and c2. */
	pos_c1 = 0;
	pos_c2 = 0;
	for (i = 0; i < size; i++) {
		/* Advance position marker until we reach a free position in c2. */
		while (pos_c2 < size && c2[pos_c2] != - 1)
			pos_c2++;
		/* c2 get the next city from p1 which is not already present. */
		if (!TestNumber(c2, size, p1[i]))
			c2[pos_c2] = p1[i];
		/* Do the same for c1. */
		while (pos_c1 < size && c1[pos_c1] != - 1)
			pos_c1++;
		/* c1 get the next city from p2 which is not already present. */
		if (!TestNumber(c1, size, p2[i]))
			c1[pos_c1] = p2[i];
	}
}	

/**
  * Crossover operator that does nothing.
  */

void fgen_crossover_noop(FgenPopulation *pop, const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	memcpy(child1, parent1, INDIVIDUAL_SIZE_IN_BYTES(pop));
	memcpy(child2, parent2, INDIVIDUAL_SIZE_IN_BYTES(pop));
}

