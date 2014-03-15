/*
    random.h -- local and inline random number generator functions.

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

// This file implements a set of inline functions used for random number generation.
// None of these are exported outside libfgen, although public random number library
// functions (implemented in random.c) make use of them. Inside libfgen, the inline
// functions can be called directly instead of exported-level functions to improve
// performance.

// This function is preferred to not be inline, because it may be called less
// frequently, and it would cause too much inlining since the RandomBit() function
// is already inlined, which covers most cases.
FGEN_LOCAL unsigned int RandomBitsNeedStorage(FgenRNG *rng, unsigned int n_bits);

FGEN_LOCAL void CalculateRandomOrder(FgenRNG *rng, int *order, int n);

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

static inline unsigned int RandomInt(FgenRNG *rng, unsigned int n) {
	return RandomIntEmpirical(rng, n);
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

static inline void RandomIntPowerOfTwoPrepareForRepeat(FgenRNG *rng, unsigned int n) {
	if (n == GetLastPowerOfTwo(rng))
		return;
	int shift = CalculateLog2(n);
	SetLastPowerOfTwoData(rng, n, shift);
}

// Repeat random integer function with the previously used power of two range.
static inline unsigned int RandomIntPowerOfTwoRepeat(FgenRNG *rng)
{
	return RandomBits(rng, GetLastPowerOfTwoShift(rng));
}

static inline void RandomIntGeneralPrepareForRepeat(FgenRNG *rng, unsigned int n) {
	if (n == GetLastGeneralRange(rng))
		return;
	int shift = CalculateBitsNeeded(n);
	SetLastGeneralRangeData(rng, n, shift);
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

