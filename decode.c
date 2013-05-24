/*
    decode.c -- bitstring decoding helper functions.

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
#include <math.h>
#include "fgen.h"
#include "parameters.h"


/**
 * Scale the 16-bit integer value stored in bitstring to a double between the given bounds.
 */

double fgen_bitstring_uint16_to_double(const unsigned char *bitstring, double domain_min, double domain_max) {
	unsigned int x_int;
	double x;
	/* Construct the 16-bit integer that the bitstring represents. */
	x_int = *(unsigned short int *)&bitstring[0];

	/* Map to the real domain [domain_min, domain_max[. */
	x = domain_min + (domain_max - domain_min) * x_int / (double)65536;
	return x;
}

/**
 * Scale the 32-bit integer value stored in bitstring to a double between the given bounds.
 */

double fgen_bitstring_uint32_to_double(const unsigned char *bitstring, double domain_min, double domain_max) {
	unsigned int x_int;
	double x;
	/* Construct the 32-bit integer that the bitstring represents. */
#if 1
	x_int = *(unsigned int *)&bitstring[0];
#else
	x_int = bitstring[0];
	x_int += (unsigned int)bitstring[1] << 8;
	x_int += (unsigned int)bitstring[2] << 16;
	x_int += (unsigned int)bitstring[3] << 24;
#endif

	/* Map to the real domain [domain_min, domain_max[. */
	x = domain_min + (domain_max - domain_min) * x_int / ((double)65536 * 65536);
	return x;
}

/**
 * Scale the 64-bit integer value stored in bitstring to a double between the given bounds.
 */

double fgen_bitstring_uint64_to_double(const unsigned char *bitstring, double domain_min, double domain_max) {
	uint64_t x_int;
	double x;
	/* Construct the 64-bit integer that the bitstring represents. */
#if 1
	x_int = *(uint64_t *)&bitstring[0];
#else
	x_int = bitstring[0];
	x_int += (uint64_t)bitstring[1] << 8;
	x_int += (uint64_t)bitstring[2] << 16;
	x_int += (uint64_t)bitstring[3] << 24;
	x_int += (uint64_t)bitstring[4] << 32;
	x_int += (uint64_t)bitstring[5] << 40;
	x_int += (uint64_t)bitstring[6] << 48;
	x_int += (uint64_t)bitstring[7] << 56;
#endif

	/* Map to the real domain [domain_min, domain_max[. */
	x = domain_min + (domain_max - domain_min) * x_int / (double)pow((double)2, 64);
	return x;
}


