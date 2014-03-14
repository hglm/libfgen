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
#include <float.h>
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

void gettimeofday(struct timeval *tv, struct timezone *tz)
{
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


/* Complementary-multiply-with-carry random number generator. */

#define PHI 0x9e3779b9

/**
 * Create a random number generator data structure and returns it. The RNG is initialized with a seed
 * of 0.
 */

FgenRNG *fgen_random_create_rng()
{
	FgenRNG *rng = (FgenRNG *)malloc(sizeof(FgenRNG));
	rng->c = 362436;
	rng->index = FGEN_RNG_STATE_SIZE - 1;
	fgen_random_seed_rng(rng, 0);
	rng->storage = 0;
	rng->storage_size = 0;
	rng->last_power_of_two = 0;
	rng->last_power_of_two_shift = 0;
	return rng;
}

/**
 * Destroy the data structure associated with an RNG.
 */

void fgen_random_destroy_rng(FgenRNG *rng)
{
	free(rng);
}

/**
 * Seed the random number generator with an unsigned integer from 0 to 2^32 - 1.
 */

void fgen_random_seed_rng(FgenRNG *rng, unsigned int seed)
{
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

unsigned int fgen_random_32(FgenRNG *rng)
{
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

// General random number generator helper functions.

static inline unsigned int GetLastPowerOfTwo(FgenRNG *rng)
{
	return rng->last_power_of_two;
}

static inline int GetLastPowerOfTwoShift(FgenRNG *rng)
{
	return rng->last_power_of_two_shift;
}

static inline unsigned int GetLastGeneralRange(FgenRNG *rng)
{
	return rng->last_general_range;
}

static inline int GetLastGeneralRangeShift(FgenRNG *rng)
{
	return rng->last_general_range_shift;
}

static inline unsigned int Random32(FgenRNG *rng)
{
	return fgen_random_32(rng);
}

static inline void SetLastPowerOfTwoData(FgenRNG *rng, unsigned int n, unsigned int shift)
{
	rng->last_power_of_two = n;
	rng->last_power_of_two_shift = shift;
}

static inline void SetLastGeneralRangeData(FgenRNG *rng, unsigned int n, unsigned int shift)
{
	rng->last_general_range = n;
	rng->last_general_range_shift = shift;
}

// Calculate floor(log2(n))). For a power of two, this is equivalent to the
// number of bits needed to represent the range 0 to n - 1. For a non-power-of-two,
// the return value is one less than the number of bits needed to represent
// the range 0 to n - 1.
static inline unsigned int CalculateLog2(unsigned int n)
{
	// Set shift to 16 if bits 15-31 are non-zero, zero otherwise.
	unsigned int shift = (((n >> 16) + 0xFFFF) & 0x10000) >> 12;
	unsigned int bits = n >> shift;
	// Add 8 to shift if bits 8-15 of the highest non-zero half-word found previously
	// are non-zero.
	unsigned char byte_shift = (((bits >> 8) + 0xFF) & 0x100) >> 5;
	shift += byte_shift;
	bits >>= byte_shift;
	// Add 4 to shift if bits 4-7 of the highest non-zero byte found previously
	// are non-zero.
	unsigned char nibble_shift = (((bits >> 4) + 0xF) & 0x10) >> 2;
	shift += nibble_shift;
	bits >>= nibble_shift;
	// Add 2 to shift if bits 2-3 of the highest non-zero nibble found previously
	// are non-zero.
	unsigned char pair_shift = (((bits >> 2) + 0x3) & 0x4) >> 1;
	shift += pair_shift;
	bits >>= pair_shift;
	// Add 1 to shift if bit 1 of the highest non-zero pair found previously
	// is non-zero.
	shift += bits >> 1;
	return shift;
}

// Calculate number of bits needed for an integer range of n (log2(n - 1) + 1).
static inline unsigned int CalculateBitsNeeded(unsigned int n)
{
	unsigned int shift = CalculateLog2(n);
	// If n is not a power of two, one more bit is needed.
	// Rely on the fact that bit 31 will be set when subtracting n from 2 ^ shift
	// and n is not power of two.
	shift += ((1 << shift) - n) >> 31;
	return shift;
}

// Calculate floor(log2(n)) when n is guaranteed to be <= 256, so that the
// return value will be <= 8. Undefined for n > 256.
static inline unsigned int CalculateLog2Max256(unsigned int n)
{
	// Set shift to 4 if bits 4-7 of the highest non-zero byte found previously
	// are non-zero.
	unsigned int shift = (((n >> 4) + 0xF) & 0x10) >> 2;
	unsigned int bits = n >> shift;
	// Add 2 to shift if bits 2-3 of the highest non-zero nibble found previously
	// are non-zero.
	unsigned char pair_shift = (((bits >> 2) + 0x3) & 0x4) >> 1;
	shift += pair_shift;
	bits >>= pair_shift;
	// Add 1 to shift if bit 1 of the highest non-zero pair found previously
	// is non-zero.
	shift += bits >> 1;
	// When n = 2^16, set shift to 16 (shift will still be zero).
	shift += (n & 0x10000) >> 12;
	return shift;
}

// Calculate number of bits needed for an integer range of n (log2(n - 1) + 1),
// when n <= 256.
static inline unsigned int CalculateBitsNeededMax256(unsigned int n)
{
	unsigned int shift = CalculateLog2Max256(n);
	// If n is not a power of two, one more bit is needed.
	// Rely on the fact that bit 31 will be set when subtracting n from 2 ^ shift
	// and n is not power of two.
	shift += ((1 << shift) - n) >> 31;
	return shift;
}

// Calculate floor(log2(n)) when n is guaranteed to be <= 2^16, so that the
// return value will be <= 16. Undefined for n > 2^16.
static inline unsigned int CalculateLog2Max65536(unsigned int n)
{
	// Set shift to 8 if bits 7-15 are non-zero.
	unsigned int shift = (((n >> 8) + 0xFF) & 0x100) >> 5;
	unsigned bits = n >> shift;
	// Add 4 to shift if bits 4-7 of the highest non-zero byte found previously
	// are non-zero.
	unsigned char nibble_shift = (((bits >> 4) + 0xF) & 0x10) >> 2;
	shift += nibble_shift;
	bits >>= nibble_shift;
	// Add 2 to shift if bits 2-3 of the highest non-zero nibble found previously
	// are non-zero.
	unsigned char pair_shift = (((bits >> 2) + 0x3) & 0x4) >> 1;
	shift += pair_shift;
	bits >>= pair_shift;
	// Add 1 to shift if bit 1 of the highest non-zero pair found previously
	// is non-zero.
	shift += bits >> 1;
	// When n = 256, set shift to 8 (shift will still be zero).
	shift += (n & 0x100) >> 5;
	return shift;
}

// Calculate number of bits needed for an integer range of n (log2(n - 1) + 1),
// when n <= 65536.
static inline unsigned int CalculateBitsNeededMax65536(unsigned int n)
{
	unsigned int shift = CalculateLog2Max65536(n);
	// If n is not a power of two, one more bit is needed.
	// Rely on the fact that bit 31 will be set when subtracting n from 2 ^ shift
	// and n is not power of two.
	shift += ((1 << shift) - n) >> 31;
	return shift;
}

int fgen_calculate_shift(unsigned int n)
{
	unsigned int shift = CalculateLog2(n);
	if (((unsigned int)1 << shift) == n)
		return shift;
	return - 1;
}

// Helper function for the inline version of RandomBits(n_bits) for when
// the storage size is known to be insufficient.

static unsigned int RandomBitsNeedStorage(FgenRNG *rng, unsigned int n_bits)
{
	unsigned int r = fgen_random_32(rng);
	// Just append the new bits at the end of the storage, possibly losing bits.
	// The higher order n_bits bits of r will be used for the return value.
	rng->storage += r << rng->storage_size;
#if FGEN_RANDOM_STORAGE_SIZE == 64
	// Nothing to do for 64-bit storage; there will always be room for at least
	// 32 bits since the maximum request size is 32.
#else
	// Adjust the storage size, limiting it to 32.
	rng->storage_size = (rng->storage_size & 32) + (((rng->storage_size & 32) >> 5) ^ 1) * rng->storage_size;
#endif
	rng->storage_size += 32 - n_bits;
	return r >> (32 - n_bits);
}

// Get n random bits. This version works for 0 <= n <= 32.
static inline unsigned int RandomBits(FgenRNG *rng, unsigned int n_bits)
{
	if (rng->storage_size < n_bits)
		return RandomBitsNeedStorage(rng, n_bits);
	unsigned int mask = ((unsigned int)1 << n_bits) - 1;
	unsigned int r = rng->storage & mask;
	rng->storage >>= n_bits;
	rng->storage_size -= n_bits;
	return r;
}

/**
 * Return a random integer value of 0 or 1.
 */

int fgen_random_2(FgenRNG *rng)
{
	return RandomBits(rng, 1);
}

/**
 * Return a random integer value from 0 to 255.
 */

int fgen_random_8(FgenRNG *rng)
{
	return RandomBits(rng, 8);
}

/**
 * Return a random integer value from 0 to 65535.
 */

int fgen_random_16(FgenRNG *rng)
{
	return RandomBits(rng, 16);
}

static inline unsigned int RandomIntEmpirical(FgenRNG *rng, unsigned int n)
{
	if (n == GetLastPowerOfTwo(rng)) {
		unsigned int shift = GetLastPowerOfTwoShift(rng);
		// Repeated bit sizes >= 20 trigger a lot of storage refills, so it is
		// faster to use Random32() and discard some bits.
		if (shift >= 20) {
			return Random32(rng) & (n - 1);
		}
		return RandomBits(rng, shift);
	}
	unsigned int shift;
	{
		shift = CalculateBitsNeeded(n);
		if (((unsigned int)1 << shift) == n) {
			SetLastPowerOfTwoData(rng, n, shift);
			return RandomBits(rng, shift);
		}
	}
	for (;;) {
		// Keep trying until the value is within the range.
		unsigned int r = RandomBits(rng, shift);
		if (r < n)
			return r;
	}
}

static inline unsigned int RandomIntEmpiricalMax256(FgenRNG *rng, unsigned int n)
{
	if (n == GetLastPowerOfTwo(rng)) {
		unsigned int shift = GetLastPowerOfTwoShift(rng);
		return RandomBits(rng, shift);
	}
	unsigned int shift;
	{
		shift = CalculateBitsNeededMax256(n);
		if (((unsigned int)1 << shift) == n) {
			SetLastPowerOfTwoData(rng, n, shift);
			return RandomBits(rng, shift);
		}
	}
	for (;;) {
		// Keep trying until the value is within the range.
		unsigned int r = RandomBits(rng, shift);
		if (r < n)
			return r;
	}
}

static inline unsigned int RandomIntEmpiricalMax65536(FgenRNG *rng, unsigned int n)
{
	if (n == GetLastPowerOfTwo(rng)) {
		unsigned int shift = GetLastPowerOfTwoShift(rng);
		return RandomBits(rng, shift);
	}
	unsigned int shift;
	{
		shift = CalculateBitsNeededMax65536(n);
		if (((unsigned int)1 << shift) == n) {
			SetLastPowerOfTwoData(rng, n, shift);
			return RandomBits(rng, shift);
		}
	}
	for (;;) {
		// Keep trying until the value is within the range.
		unsigned int r = RandomBits(rng, shift);
		if (r < n)
			return r;
	}
}

/* The following functions use the RNG implementation but do not not differ between implementations. */

/**
 * Randomize the seed of the random number generator with a value from the system timer.
 */

void fgen_random_seed_with_timer(FgenRNG *rng)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	/* The multiplication by 1000000 will overflow, but that is not important. */
	fgen_random_seed_rng(rng, tv.tv_sec * 1000000 + tv.tv_usec);
}

/**
 * Return an integer from 0 to n - 1.
 */

int fgen_random_n(FgenRNG *rng, int n)
{
	return RandomIntEmpirical(rng, n);
}

unsigned int fgen_random_n_max_65536(FgenRNG *rng, unsigned int n)
{
	return RandomIntEmpiricalMax65536(rng, n);
}

unsigned int fgen_random_n_max_256(FgenRNG *rng, unsigned int n)
{
	return RandomIntEmpiricalMax256(rng, n);
}

// Any power of two range 1 <= n <= 256.
static inline unsigned int RandomIntPowerOfTwoMax256(FgenRNG *rng, unsigned int n)
{
	// Use the calculation method.
	unsigned int shift = CalculateLog2Max256(n);
	return RandomBits(rng, shift);
}

// Any power of two range 1 <= n <= 2^16.
static inline unsigned int RandomIntPowerOfTwoMax65536(FgenRNG *rng, unsigned int n)
{
	// Use the calculation method.
	unsigned int shift = CalculateLog2Max65536(n);
	return RandomBits(rng, shift);
}

// Random integer for general power-of-two range from 1 up to (1 << 30).
static inline unsigned int RandomIntPowerOfTwo(FgenRNG *rng, unsigned int n)
{
	if (n == GetLastPowerOfTwo(rng)) {
		return RandomBits(rng, GetLastPowerOfTwoShift(rng));
	}
	unsigned int shift;
	// Use the calculation method.
	shift = CalculateLog2(n);
	SetLastPowerOfTwoData(rng, n, shift);
	return RandomBits(rng, shift);
}

// Return random integer in any power of two range with specified shift
// (log2(n)). Different from RandomBits() because it caches the power of
// two range value.
static inline unsigned int RandomIntPowerOfTwoWithShift(FgenRNG *rng, unsigned int shift)
{
	SetLastPowerOfTwoData(rng, 1 << shift, shift);
	return RandomBits(rng, shift);
}

// Repeat random integer function with the previously used power of two range.
static inline unsigned int RandomIntPowerOfTwoRepeat(FgenRNG *rng)
{
	return RandomBits(rng, GetLastPowerOfTwoShift(rng));
}

// Repeat random integer function with the non-power-of-two range previously set
// with RandomIntGeneralPrepareForRepeat().
static inline unsigned int RandomIntGeneralRepeat(FgenRNG *rng)
{
	for (;;) {
		// Keep trying until the value is within the range.
		unsigned int r = RandomBits(rng, GetLastGeneralRangeShift(rng));
		if (r < GetLastGeneralRange(rng))
			return r;
	}
}

// Additional fgen_random() functions for powers of two (to be implemented).

unsigned int fgen_random_n_power_of_two(FgenRNG *rng, unsigned int n)
{
	return RandomIntPowerOfTwo(rng, n);
}

unsigned int fgen_random_n_power_of_two_max_65536(FgenRNG *rng, unsigned int n)
{
	return RandomIntPowerOfTwoMax65536(rng, n);
}

unsigned int fgen_random_n_power_of_two_max_256(FgenRNG *rng, unsigned int n)
{
	return RandomIntPowerOfTwoMax256(rng, n);
}

unsigned int fgen_random_n_power_of_two_with_shift(FgenRNG *rng, unsigned int shift)
{
	return RandomIntPowerOfTwoWithShift(rng, shift);
}

unsigned int fgen_random_n_power_of_two_repeat(FgenRNG *rng)
{
	return RandomIntPowerOfTwoRepeat(rng);
}

void fgen_random_n_prepare_for_power_of_two(FgenRNG *rng, unsigned int n)
{
	if (n == GetLastPowerOfTwo(rng))
		return;
	int shift = CalculateLog2(n);
	SetLastPowerOfTwoData(rng, n, shift);
}

unsigned int fgen_random_n_general_repeat(FgenRNG *rng)
{
	return RandomIntGeneralRepeat(rng);
}

void fgen_random_n_general_prepare_for_repeat(FgenRNG *rng, unsigned int n)
{
	if (n == GetLastGeneralRange(rng))
		return;
	int shift = CalculateBitsNeeded(n);
	SetLastGeneralRangeData(rng, n, shift);
}

/**
 * Return a random double from 0 to range (exclusive).
 */

double fgen_random_d(FgenRNG *rng, double range)
{
	// Return high-precision double.
	// Scaling the 2^32 integers to [0, 1) (which has good precision in the
	// double format) maintains precision when adding up the lower and higher
	// order components.
	return ((double)Random32(rng) * (1.0d / pow(2.0d, 32)) +
	        (double)Random32(rng) * (1.0d / pow(2.0d, 64))) * range;
}

/**
 * Return a random float from 0 to range (exclusive).
 */

float fgen_random_f(FgenRNG *rng, float range)
{
	// Return high-precision random float.
	return (float)((uint64_t)Random32(rng) << 32) * (1.0f / powf(2.0f, 64)) * range;
}

float fgen_random_f_low_precision(FgenRNG *rng, float range)
{
	return (float)Random32(rng) * (1.0f / powf(2.0f, 32)) * range;
}

float fgen_random_f_very_low_precision(FgenRNG *rng, float range)
{
	return (float)(unsigned short)RandomBits(rng, 16) * (1.0f / powf(2.0f, 16)) * range;
}

float fgen_random_d_low_precision(FgenRNG *rng, float range)
{
	return (double)Random32(rng) * (1.0f / pow(2.0d, 32)) * range;
}

float fgen_random_d_high_precision(FgenRNG *rng, float range)
{
	if (range <= 1.00001d)
		// When range <= 1.0, the standard high-precision function is already
		// optimal.
		return fgen_random_d(rng, range);
	const double high_value = DBL_MAX;
	// Use the identity exp(x + y) = exp(x) * exp(y).
	// Scale the 32 bit random integers r0 and r1 so that
	// 0 <= exp(scale * r0 * exp(scale * r1 * 2^32) <= high_value,
	// where high_value is near the greatest representable double value.
	// That is, exp(scale * r0 + scale * r1 * 2^32) <= high_value,
	// <-> scale * r0 + scale * r1 * 2^32 <= log(high_value)
	// so that scale = log(high_value) / pow(2.0d, 32).
	const double scale_factor0 = log(high_value) / pow(2.0d, 32);
	const double scale_factor1 = log(high_value) / pow(2.0d, 64);
	return log(
	               exp((double)Random32(rng) * scale_factor0) *
	               exp((double)Random32(rng) * scale_factor1)
	       ) * (range / log(high_value));
}

/**
 * Return a random double from min_bound to max_bound (exclusive).
 */

double fgen_random_from_range_d(FgenRNG *rng, double min_bound, double max_bound)
{
	return min_bound + fgen_random_d(rng, max_bound - min_bound);
}

/* Calculate a random permutation of the numbers 0 to (n - 1). */

void CalculateRandomOrder(FgenRNG *rng, int *order, int n)
{
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


