/*
    random.c -- random number generator functions.

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
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifdef __GNUC__
#include <sys/time.h>
#include <pthread.h>
#else
#include <time.h>
#include <windows.h>
#endif
#include "fgen.h"
#include "random.h"
#include "error.h"

#ifndef __GNUC__

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

void gettimeofday(struct timeval *tv, struct timezone *tz) {
	FILETIME ft;
	uint64_t tmpres = 0;
	if (NULL != tv) {
		GetSystemTimeAsFileTime(&ft);
		tmpres |= ft.dwHighDateTime;
		tmpres <<= 32;
		tmpres |= ft.dwLowDateTime;
		tmpres /= 10;  /*convert into microseconds*/
		/*converting file time to unix epoch*/
		tmpres -= DELTA_EPOCH_IN_MICROSECS; 
		tv->tv_sec = (long)(tmpres / 1000000UL);
		tv->tv_usec = (long)(tmpres % 1000000UL);
	}
 }

#endif

/* The random number generator to use. glibc slows down with mutex blocking when running concurrent threads that
 * call random() */

// #define RANDOM_USE_GLIBC_ONLY
// #define RANDOM_USE_GLIBC_OPTIMIZED
#define RANDOM_USE_CMWC

#ifdef RANDOM_USE_CMWC

/* Complementary-multiply-with-carry random number generator. */

#define PHI 0x9e3779b9

/**
 * Create a random number generator data structure and returns it. The RNG is initialized with a seed
 * of 0.
 */

FgenRNG *fgen_random_create_rng() {
	FgenRNG *rng = (FgenRNG *)malloc(sizeof(FgenRNG));
	rng->c = 362436;
	rng->index = FGEN_RNG_STATE_SIZE - 1;
	fgen_random_seed_rng(rng, 0);
	rng->storage = 0;
	rng->storage_size = 0;
	rng->last_random_n_power_of_2 = - 1;
	return rng;
}

/**
 * Destroy the data structure associated with an RNG.
 */

void fgen_random_destroy_rng(FgenRNG *rng) {
	free(rng);
}

/**
 * Seed the random number generator with an unsigned integer from 0 to 2^32 - 1.
 */

void fgen_random_seed_rng(FgenRNG *rng, unsigned int seed) {
        int i;
        rng->Q[0] = seed;
        rng->Q[1] = seed + PHI;
        rng->Q[2] = seed + PHI + PHI;
        for (i = 3; i < FGEN_RNG_STATE_SIZE; i++)
                rng->Q[i] = rng->Q[i - 3] ^ rng->Q[i - 2] ^ PHI ^ i;
}
 
/**
 * Return a random integer value from 0 to 2^32 - 1;
 */

unsigned int fgen_random_32(FgenRNG *rng) {
        uint64_t t, a = 18782LL;
        unsigned int x, r = 0xfffffffe;
        rng->index = (rng->index + 1) & (FGEN_RNG_STATE_SIZE - 1);
        t = a * rng->Q[rng->index] + rng->c;
        rng->c = (t >> 32);
        x = t + rng->c;
        if (x < rng->c) {
                x++;
                rng->c++;
        }
        return (rng->Q[rng->index] = r - x);
}

/**
 * Return a random integer value of 0 or 1.
 */

int fgen_random_2(FgenRNG *rng) {
	unsigned int r;
	if (rng->storage_size > 0) {
		int bit;
		bit = rng->storage & 0x1;
		rng->storage >>= 1;
		rng->storage_size--;
		return bit;
	}
	r = fgen_random_32(rng);
	rng->storage = (r & 0xFFFFFFFE) >> 1;
	rng->storage_size = 31;
	return r & 0x1;
}

/**
 * Return a random integer value from 0 to 255.
 */

int fgen_random_8(FgenRNG *rng) {
	unsigned int r;
	if (rng->storage_size >= 8) {
		r = rng->storage & 0xFF;
		rng->storage >>= 8;
		rng->storage_size -= 8;
		return r;
	}
	r = fgen_random_32(rng);
	rng->storage += ((r & 0xFFFFFF00) >> 8) << rng->storage_size;
	rng->storage_size += 24;
	return r & 0xFF;
}

/**
 * Return a random integer value from 0 to 65535.
 */

int fgen_random_16(FgenRNG *rng) {
	unsigned int r;
	if (rng->storage_size >= 16) {
		r = rng->storage & 0xFFFF;
		rng->storage >>= 16;
		rng->storage_size -= 16;
		return r;
	}
	r = fgen_random_32(rng);
	rng->storage += ((r & 0xFFFF0000) >> 16) << rng->storage_size;
	rng->storage_size += 16;
	return r & 0xFFFF;
}

#endif


#ifdef RANDOM_USE_GLIBC_ONLY

/**
 * Create a random number generator data structure and returns it. The RNG is initialized with a seed
 * of 0.
 */

FgenRNG *fgen_random_create_rng() {
	FgenRNG *rng = malloc(sizeof(FgenRNG));
	fgen_random_seed_rng(rng, 0);
//	rng->storage = 0;
//	rng->storage_size = 0;
	return rng;
}

/**
 * Destroy the data structure associated with an RNG.
 */

void fgen_random_destroy_rng(FgenRNG *rng) {
	free(rng);
}

/**
 * Seed the random number generator with an unsigned integer from 0 to 2^32 - 1.
 */

void fgen_random_seed_rng(FgenRNG *rng, unsigned int seed) {
	srandom(seed);
}
 
/**
 * Return a random integer value from 0 to 2^32 - 1;
 */

unsigned int fgen_random_32(FgenRNG *rng) {
	unsigned int r = (unsigned int)random() | (((unsigned int)random() & 0x1) << 31);
	return r;
}

/**
 * Return a random integer value of 0 or 1.
 */

int fgen_random_2(FgenRNG *rng) {
	int r = random() & 1;
	return r;
}

/**
 * Return a random integer value from 0 to 255.
 */

int fgen_random_8(FgenRNG *rng) {
	int r = random() & 0xFF;
	return r;
}

/**
 * Return a random integer value from 0 to 65535.
 */

int fgen_random_16(FgenRNG *rng) {
	int r = random() & 0xFFFF;
	return r;
}

#endif

#ifdef RANDOM_USE_GLIBC_OPTIMIZED

/*
 * The random() function usually returns numbers between 0 and MAXINT31;
 * to get a 32-bit number, two calls are needed.
 *
 * We try to be smart and keep a collection of left-over 'random bits' to
 * reduce the number of calls to random() (for the case RAND_MAX == MAXINT31,
 * which is most common).
 *
 */

#define MAXINT31 ((unsigned)(1 << 31) - 1)

/**
 * Create a random number generator data structure and returns it. The RNG is initialized with a seed
 * of 0.
 */

FgenRNG *fgen_random_create_rng() {
	FgenRNG *rng = malloc(sizeof(FgenRNG));
	fgen_random_seed_rng(rng, 0);
	rng->storage = 0;
	rng->storage_size = 0;
	rng->last_random_n_power_of_2 = - 1;
	return rng;
}

/**
 * Destroy the data structure associated with an RNG.
 */

void fgen_random_destroy_rng(FgenRNG *rng) {
	free(rng);
}

/**
 * Seed the random number generator with an unsigned integer from 0 to 2^32 - 1.
 */

void fgen_random_seed_rng(FgenRNG *rng, unsigned int seed) {
	srandom(seed);
}
 
/**
 * Return a random integer value from 0 to 2^32 - 1;
 */

unsigned int fgen_random_32(FgenRNG *rng) {
	unsigned int r1, r2;
	if (rng->storage_size > 0) {
		/* Need one storage bit. */
		unsigned int bit;
		bit = (rng->storage & 0x1) << 31;
		rng->storage >>= 1;
		rng->storage_size--;
		return random() + bit;
	}
	/*
	 * No storage bits; do two random()'s and put the 15 + 15 = 30 extra
	 * bits in storage.
	 */
	r1 = random();
	rng->storage = (unsigned int)(r1 & 0x7FFF0000) >> 16;
	r2 = random();
	rng->storage += (unsigned int)(r2 & 0x7FFF0000) >> 1;
	rng->storage_size = 30;
	return (r1 & 0xFFFF) + ((unsigned int)(r2 & 0xFFFF) << 16);
}

/**
 * Return a random integer value of 0 or 1.
 */

int fgen_random_2(FgenRNG *rng) {
	unsigned int r;
	if (rng->storage_size > 0) {
		int bit;
		bit = rng->storage & 0x1;
		rng->storage >>= 1;
		rng->storage_size--;
		return bit;
	}
	r = fgen_random_32(rng);
	rng->storage = (r & 0xFFFFFFFE) >> 1;
	rng->storage_size = 31;
	return r & 0x1;
}

/**
 * Return a random integer value from 0 to 255.
 */

int fgen_random_8(FgenRNG *rng) {
	unsigned int r;
	if (rng->storage_size >= 8) {
		r = rng->storage & 0xFF;
		rng->storage >>= 8;
		rng->storage_size -= 8;
		return r;
	}
	r = fgen_random_32(rng);
	rng->storage += ((r & 0xFFFFFF00) >> 8) << rng->storage_size;
	rng->storage_size += 24;
	return r & 0xFF;
}

/**
 * Return a random integer value from 0 to 65535.
 */

int fgen_random_16(FgenRNG *rng) {
	unsigned int r;
	if (rng->storage_size >= 16) {
		r = rng->storage & 0xFFFF;
		rng->storage >>= 16;
		rng->storage_size -= 16;
		return r;
	}
	r = fgen_random_32(rng);
	rng->storage += ((r & 0xFFFF0000) >> 16) << rng->storage_size;
	rng->storage_size += 16;
	return r & 0xFFFF;
}

#endif

/* The following functions use the RNG implementation but do not not differ between implementations. */

/**
 * Randomize the seed of the random number generator with a value from the system timer.
 */

void fgen_random_seed_with_timer(FgenRNG *rng) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	/* The multiplication by 1000000 will overflow, but that is not important. */
	fgen_random_seed_rng(rng, tv.tv_sec * 1000000 + tv.tv_usec);
}

static unsigned char power_of_2_table[257] = {
	255, 0, 1, 255, 2, 255, 255, 255, 3, 255, 255, 255, 255 ,255, 255, 255, 
	4, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	5, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	6, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255, 
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 ,255, 255, 255,
	8
};

static int fgen_random_get_bits(FgenRNG *rng, int n) {
	unsigned int r;
	unsigned int mask = (1 << n) - 1;
	if (rng->storage_size >= n) {
		r = rng->storage & mask;
		rng->storage >>= n;
		rng->storage_size -= n;
		return r;
	}
	r = fgen_random_32(rng);
	rng->storage |= ((r & (0xFFFFFFFF ^ mask)) >> n) << rng->storage_size;
	rng->storage_size += 32 - n;
	return r & mask;
}

/**
 * Return an integer from 0 to n - 1.
 */

int fgen_random_n(FgenRNG *rng, int n) {
	// Fast path to most common occurence of repeating power of two.
	if (n == rng->last_random_n_power_of_2) {
		return fgen_random_get_bits(rng, rng->last_random_n_power_of_2_bit_count);
	}
	// Optimize power of two.
	if (n <= 256) {
		int bit_count = power_of_2_table[n];
		if (bit_count != 255) {
			rng->last_random_n_power_of_2 = n;
			rng->last_random_n_power_of_2_bit_count = bit_count;
			return fgen_random_get_bits(rng, bit_count);
		}
		return fgen_random_16(rng) % n;
	}
	if (n <= 65536 && (n & 0xFF) == 0) {
		int bit_count = power_of_2_table[n >> 8];
		if (bit_count != 255) {
			rng->last_random_n_power_of_2 = n;
			rng->last_random_n_power_of_2_bit_count = bit_count + 8;
			return fgen_random_get_bits(rng, bit_count + 8);
		}
	}
	return fgen_random_32(rng) % n;
}

/** 
 * Return a random double from 0 to range (exclusive).
 */

double fgen_random_d(FgenRNG *rng, double range) {
    return (double)fgen_random_32(rng) * range / ((uint64_t)1 << 32);
}

/**
 * Return a random float from 0 to range (exclusive).
 */

float fgen_random_f(FgenRNG *rng, float range) {
    return (float)fgen_random_32(rng) * range / ((uint64_t)1 << 32);
}

/**
 * Return a random double from min_bound to max_bound (exclusive).
 */

double fgen_random_from_range_d(FgenRNG *rng, double min_bound, double max_bound) {
	return min_bound + fgen_random_d(rng, max_bound - min_bound);
}

/* Calculate a random permutation of the numbers 0 to (n - 1). */

void CalculateRandomOrder(FgenRNG *rng, int *order, int n) {
	int i;
	for (i = 0; i < n; i++)
		order[i] = i;
	for (i = 0; i < n; i++) {
		int j, t;
		/* Swap element i with random element j. */
		j = fgen_random_n(rng, n);
		t = order[i];
		order[i] = order[j];
		order[j] = t;
	}
}


