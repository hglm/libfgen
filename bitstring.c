/*
    bitstring.c -- functions for handling bitstrings.

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
#include <stdint.h>
#include <limits.h>
#include <fgen.h>
#include <random.h>

#if __WORDSIZE == 64
#define POINTER_TO_UINT(x) ((unsigned long int)(x))
#else
#define POINTER_TO_UINT(x) ((unsigned int)(x))
#endif

void CreateRandomBitstring(FgenRNG *rng, unsigned char *bitstring, int size_in_bytes) {
	while (size_in_bytes >= 1 && (POINTER_TO_UINT(bitstring) & 3) != 0) {
		*bitstring = RandomBits(rng, 8);
		size_in_bytes--;
		bitstring++;
	}
	while (size_in_bytes >= 4) {
		*(unsigned int *)bitstring = fgen_random_32(rng);
		size_in_bytes -= 4;
		bitstring += 4;
	}
	while (size_in_bytes >= 1) {
		*bitstring = RandomBits(rng, 8);
		size_in_bytes--;
		bitstring++;
	}
}

/**
 * Set the bit at position offset in bitstring to a random value (1 or 0). offset must be >= 0 but may
 * be greater than 7.
 */

static void fgen_set_random_bit(FgenRNG *rng, unsigned char *bitstring, const int offset) {
    bitstring[offset / 8] &= 0xFF ^ (1 << (offset & 7));
    bitstring[offset / 8] |= fgen_random_2(rng) << (offset & 7);
}

/**
 * Set a number of bits at position offset in bitstring each to a random value (1 or 0). offset must be >= 0
 * but may be greater than 7.
 */

void fgen_set_random_bitstring(FgenRNG *rng, unsigned char *bitstring, const int offset, const int nu_bits) {
    int i;
    if (nu_bits == 1) {
         fgen_set_random_bit(rng, bitstring, offset);
         return;
    }
    /* Fast case when offset and nu_bits are a multiple of 8. */
    if ((offset & 7) == 0 && (nu_bits & 7) == 0) {
        CreateRandomBitstring(rng, bitstring + offset / 8, nu_bits / 8);
        return;
    }
    for (i = offset; i < offset + nu_bits; i++)
        fgen_set_random_bit(rng, bitstring, i);
}

/**
 * Mutate (toggle) the bit at offset bitnumber in bitstring. bitnumber must be >= 0 but may be greater than 7.
 */

void fgen_mutate_bit(unsigned char *bitstring, int bitnumber) {
	bitstring[bitnumber / 8] ^= 1 << (bitnumber & 7);
}

/**
 * Returns the value (1 or 0) of the bit at offset bitnumber in bitstring. bitnumber must be >= 0.
 */

int fgen_get_bit(const unsigned char *bitstring, int bitnumber) {
	return (bitstring[bitnumber / 8] & (1 << (bitnumber & 7))) >> (bitnumber & 7);
}


/**
 * Copy the part of the source bitstring at offset bitoffset of length number_of_bits to the same position in the
 * destination bitstring.
 */

void fgen_copy_partial_bitstring(const unsigned char *src, unsigned char *dest, int bitoffset, int number_of_bits) {
	int i, offset;
	if (number_of_bits == 0)
		return;
	offset = bitoffset / 8;
	if ((bitoffset & 7) != 0) {
		for (i = bitoffset & 7; i < 8; i++) {
			dest[offset] &= ~(1 << i);
			dest[offset] |= src[offset] & (1 << i);
			number_of_bits--;
			if (number_of_bits == 0)
				return;
		}
		offset++;
	}
	while (number_of_bits >= 8 && (POINTER_TO_UINT(&dest[offset]) & 3) != 0) {
		dest[offset] = src[offset];
		number_of_bits -= 8;
		offset++;
	}
	while (number_of_bits >= 32) {
		*(unsigned int *)&dest[offset] = *(unsigned int *)&src[offset];
		number_of_bits -= 32;
		offset += 4;
	}
	while (number_of_bits >= 8) {
		dest[offset] = src[offset];
		number_of_bits -= 8;
		offset++;
	}
	for (i = 0; i < number_of_bits; i++) {
		dest[offset] &= ~(1 << i);
		dest[offset] |= src[offset] & (1 << i);
	}
}

