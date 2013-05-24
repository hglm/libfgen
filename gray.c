/*
    gray.c -- Gray coding functions.

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

#include <stdint.h>
#include "fgen.h"
#include "parameters.h"

/* Functions for Gray coding conversion. */
/* Inspired by Usenet comp.ai.genetic FAQ. */


static void ConvertBinaryToGray(const unsigned char *src, unsigned char *dest, int size_in_bytes) {
	int bitoffset;
	unsigned char *srcp, *destp;
	unsigned char savedbit;
	bitoffset = size_in_bytes * 8 - 1;
	srcp = (unsigned char *)src + bitoffset / 8;
	destp = (unsigned char *)src + bitoffset / 8;
	savedbit = *srcp & (1 << 7);
	*destp &= ~(1 << 7);
	*destp |= savedbit;
	savedbit >>= 7;
	bitoffset--;
	for (;;) {
		int i;
		if (bitoffset < 0)
			break;
		i = bitoffset & 7;
		/*
		 * For second and following bytes, to derive the
		 * high order bit for the gray code, the
		 * saved last bit in the binary string from the
		 * previous byte is used.
		 */
		while (i >= 0) {
			unsigned char bit;
			bit = *srcp & (1 << i);
			*destp &= ~(1 << i);
			*destp |= (savedbit << i) ^ bit;
			savedbit = bit >> i;
			i--;
		}
		bitoffset = (bitoffset & (~7)) - 1;
		srcp--;
		destp--;
	}
}

// General Gray-to-binary decoding function that is very slow.

#define SetBit(str, offset, bit) \
	{ \
		unsigned char mask; \
		unsigned char shift; \
		shift = (offset) & 7; \
		mask = 1 << shift; \
		str[(offset) / 8] &= 0xFF - mask; \
		str[(offset) / 8] |= (bit) << shift; \
	}

#define GetBit(str, offset) \
	((str[(offset) / 8] & (1 << ((offset) & 7))) >> ((offset) & 7))

static void ConvertGrayToBinary(const unsigned char *src, unsigned char *dest, int size_in_bytes) {
	int i;
	i = size_in_bytes * 8 - 1;
	SetBit(dest, i, GetBit(src, i));
	i--;
	while (i >= 0) {
		SetBit(dest, i, GetBit(dest, i + 1) ^ GetBit(src, i));
		i--;
	}
#if 0
	int bitoffset;
	unsigned char *srcp, *destp;
	unsigned char savedbit;
	bitoffset = size_in_bytes * 8 - 1;
	srcp = src + bitoffset / 8;
	destp = src + bitoffset / 8;
	savedbit = *srcp & (1 << 7);
	*destp &= ~(1 << 7);
	*destp |= savedbit;
	savedbit >>= 7;
	bitoffset--;
	for (;;) {
		int i;
		if (bitoffset < 0)
			break;
		i = bitoffset & 7;
		/*
		 * For second and following bytes, to derive the
		 * high order bit for the binary code, the
		 * saved last bit in the binary string from the
		 * previous byte is used.
		 */
		while (i >= 0) {
			unsigned char bit;
			bit = (*srcp & (1 << i)) ^ (savedbit << i);
			*destp &= ~(1 << i);
			*destp |= bit;
			savedbit = bit >> i;
			i--;
		}
		bitoffset = (bitoffset & (~7)) - 1;
		srcp--;
		destp--;
	}
#endif
}

static unsigned int ConvertGrayToBinary8(unsigned int num) {
	unsigned int numBits = 8;
	unsigned int shift;
	for (shift = 1; shift < numBits; shift = 2 * shift)
	        num = num ^ (num >> shift);
	return num;
}

static unsigned int ConvertGrayToBinary16(unsigned int num) {
	unsigned int numBits = 16;
	unsigned int shift;
	for (shift = 1; shift < numBits; shift = 2 * shift)
	        num = num ^ (num >> shift);
	return num;
}

static unsigned int ConvertGrayToBinary32(unsigned int num) {
	unsigned int numBits = 8 * sizeof(num);
	unsigned int shift;
	for (shift = 1; shift < numBits; shift = 2 * shift)
	        num = num ^ (num >> shift);
	return num;
}

static uint64_t ConvertGrayToBinary64(uint64_t num) {
	unsigned int numBits = 8 * sizeof(num);
	unsigned int shift;
	for (shift = 1; shift < numBits; shift = 2 * shift)
	        num = num ^ (num >> shift);
	return num;
}

// Optimized decoding functions for common data element sizes.

static void decode_from_gray_8(const FgenPopulation *pop, const unsigned char *src, unsigned char *dest) {
	int n = INDIVIDUAL_SIZE_IN_BYTES(pop);
	for (int i = 0; i < n; i++)
		dest[i] = ConvertGrayToBinary8(src[i]);
}

static void decode_from_gray_16(const FgenPopulation *pop, const unsigned short *src, unsigned short *dest) {
	int n = INDIVIDUAL_SIZE_IN_BYTES(pop) / 2;
	for (int i = 0; i < n; i++)
		dest[i] = ConvertGrayToBinary16(src[i]);
}

static void decode_from_gray_32(const FgenPopulation *pop, const unsigned int *src, unsigned int *dest) {
	int n = INDIVIDUAL_SIZE_IN_BYTES(pop) / 4;
	for (int i = 0; i < n; i++)
		dest[i] = ConvertGrayToBinary32(src[i]);
}

static void decode_from_gray_64(const FgenPopulation *pop, const uint64_t *src, uint64_t *dest) {
	int n = INDIVIDUAL_SIZE_IN_BYTES(pop) / 8;
	for (int i = 0; i < n; i++)
		dest[i] = ConvertGrayToBinary64(src[i]);
}

/**
 * Decode the source bitstring from Gray-code and store it in the destination bitstring.
 */

void fgen_decode_from_gray(const FgenPopulation *pop, const unsigned char *src_bitstring,
unsigned char *dest_bitstring) {
	if (pop->data_element_size == 1 || pop->data_element_size == pop->individual_size_in_bits)
		ConvertGrayToBinary(src_bitstring, dest_bitstring, INDIVIDUAL_SIZE_IN_BYTES(pop));
	else {
		switch (pop->data_element_size) {
		case 8 :
			decode_from_gray_8(pop, src_bitstring, dest_bitstring);
			return;
		case 16 :
			decode_from_gray_16(pop, (const unsigned short *)src_bitstring, (unsigned short *)dest_bitstring);
			return;
		case 32 :
			decode_from_gray_32(pop, (const unsigned int *)src_bitstring, (unsigned int *)dest_bitstring);
			return;
		case 64 :
			decode_from_gray_64(pop, (const uint64_t *)src_bitstring, (uint64_t *)dest_bitstring);
			return;
		default :
			/* For each data element */
			for (int i = 0; i < INDIVIDUAL_SIZE_IN_BYTES(pop) * 8 / pop->data_element_size; i++)
				ConvertGrayToBinary(src_bitstring + i * pop->data_element_size / 8, dest_bitstring +
					i * pop->data_element_size / 8, pop->data_element_size / 8);
			break;
		}
	}
}

/**
 * Encode the source bitstring to Gray-code and store it in the destination bitstring.
 */

void fgen_encode_to_gray(const FgenPopulation *pop, const unsigned char *src_bitstring,
unsigned char *dest_bitstring) {
	if (pop->data_element_size == 1 || pop->data_element_size == pop->individual_size_in_bits)
		ConvertBinaryToGray(src_bitstring, dest_bitstring, INDIVIDUAL_SIZE_IN_BYTES(pop));
	else {
		int i;
		/* For each data element */
		for (i = 0; i < INDIVIDUAL_SIZE_IN_BYTES(pop) * 8 / pop->data_element_size; i++)
			ConvertBinaryToGray(src_bitstring + i * pop->data_element_size / 8, dest_bitstring +
				i * pop->data_element_size / 8, pop->data_element_size / 8);
	}
}

