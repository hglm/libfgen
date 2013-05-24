/*
    ga.c -- genetic algorithm/programming module of gp, a linear genetic programming example.

    fgen -- Library for optimization using a genetic algorithm or particle swarm optimization.
    Copyright 2013, Harm Hanemaaijer

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
#include <limits.h>
#include <float.h>
#ifdef __GNUC__
#include <fenv.h>
#endif
#include <errno.h>
#include <fgen.h>
#include "gp.h"

#ifndef __GNUC__
#define INFINITY DBL_MAX
#endif

#ifdef _WIN32
#define random() rand()
#endif

#if defined(_WIN32) && !defined(__GNUC__)

#define isnan(x) _isnan(x)

static int custom_round(float v, int direction) {
	switch (direction) {
	case 0 :
		if (v > 0)	
			return floor(v);
		else
			return ceil(v);
	case 1 :
		return floor(v + 0.5);
	case 2 :
		return ceil(v);
	case 3 :
		return floor(v);
	}
}

#else

static int custom_round(float v, int direction) {
	switch (direction) {
	case 0 :
		fesetround(FE_TOWARDZERO);
		break;
	case 1 :
		fesetround(FE_TONEAREST);
		break;
	case 2 :
		fesetround(FE_UPWARD);
		break;
	case 3 :
		fesetround(FE_DOWNWARD);
		break;
	}
	return nearbyintf(v);
}

#endif

static char *register_str[32] = { "r0", "r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8", "r9", "r10",
	"r11", "r12", "r13", "r14", "r15", "r16", "r17", "r18", "r19", "r20",
	"r21", "r22", "r23", "r24", "r25", "r26", "r27", "r28", "r29", "r30", "r31" };

static char *fregister_str[32] = { "f0", "f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "f10",
	"f11", "f12", "f13", "f14", "f15", "f16", "f17", "f18", "f19", "f20",
	"f21", "f22", "f23", "f24", "f25", "f26", "f27", "f28", "f29", "f30", "f31" };

static void print_program(unsigned int *program) {
	int ip = 0;
	int nu_instructions = (int)program[max_instructions];
	int operator_and = 31;
	if (only_integer_math)
		operator_and = 15;
	for (;;) {
		int src_register1 = (program[ip] >> 17) & (nu_registers - 1);
		int src_register2 = (program[ip] >> 22) & (nu_registers - 1);
		int dest_register = (program[ip] >> 27) & (nu_registers - 1);
		int opcode = program[ip] & operator_and;
		if (only_floating_point_math && ((opcode >= 0 && opcode <= 8) || (opcode >= 11 && opcode <= 14)))
			opcode |= 16;
		switch (opcode) {
		case 0 :	// mov
			printf("	mov %s, %s\n", register_str[src_register1], register_str[dest_register]);
			break;
		case 1 :	// add
			printf("	add %s, %s, %s\n", register_str[src_register1], register_str[src_register2],
				register_str[dest_register]);
			break;
		case 2 :	// sub
			printf("	sub %s, %s, %s\n", register_str[src_register1], register_str[src_register2],
				register_str[dest_register]);
			break;
		case 3 :	// mul
			printf("	mul %s, %s, %s\n", register_str[src_register1], register_str[src_register2],
				register_str[dest_register]);
			break;
		case 4 :	// div
			printf("	div %s, %s, %s\n", register_str[src_register1], register_str[src_register2],
				register_str[dest_register]);
			break;
		case 5 :	// and
			printf("	and %s, %s, %s\n", register_str[src_register1], register_str[src_register2],
				register_str[dest_register]);
			break;
		case 6 :	// or
			printf("	or %s, %s, %s\n", register_str[src_register1], register_str[src_register2],
				register_str[dest_register]);
			break;
		case 7 :	// seteq
			printf("	seteq %s, %s, %s\n", register_str[src_register1], register_str[src_register2],
				register_str[dest_register]);
			break;
		case 8 :	// setgr
			printf("	setgr %s, %s, %s\n", register_str[src_register1], register_str[src_register2],
				register_str[dest_register]);
			break;
		case 9 :	// jmpz
			printf("	jmpz %s, offset(%d)\n", register_str[dest_register], (program[ip] >> 5) & 31);
			break;
		case 10: 	// jmpnz
			printf("	jmpnz %s, offset(%d)\n", register_str[dest_register], (program[ip] >> 5) & 31);
			break;
		case 11 : 	// movsew
			if (program[ip] & (1 << 20))	// Sign bit is set.
				printf("	movsew $-%d, %s\n", ((program[ip] >> 5) & 0x7FFF),
					register_str[dest_register]);
			else
				printf("	movsew $%d, %s\n", ((program[ip] >> 5) & 0xFFFF),
					register_str[dest_register]);
			break;
		case 12 : 	// movlo
			printf("	movlo $0x%04x, %s\n", (program[ip] >> 5) & 0xFFFF, register_str[dest_register]);
			break;
		case 13 : 	// movhi
			printf("	movhi $0x%04x, %s\n", (program[ip] >> 5) & 0xFFFF, register_str[dest_register]);
			break;
		case 14 :	// movseb
			if (program[ip] & (1 << 12))	// Sign bit is set.
				printf("	movseb $-%d, %s\n", 0x80 - ((program[ip] >> 5) & 0x7F),
					register_str[dest_register]);
			else
				printf("	movseb $%d, %s\n", ((program[ip] >> 5) & 0x7F), register_str[dest_register]);
			break;
		case 15 :	// end
			printf("	end\n");
			return;
		case 16 :	// movf
			printf("	movf %s, %s\n", fregister_str[src_register1], fregister_str[dest_register]);
			break;
		case 17 :	// addf
			printf("	addf %s, %s, %s\n", fregister_str[src_register1], fregister_str[src_register2],
				fregister_str[dest_register]);
			break;
		case 18 :	// subf
			printf("	subf %s, %s, %s\n", fregister_str[src_register1], fregister_str[src_register2],
				fregister_str[dest_register]);
			break;
		case 19 :	// mulf
			printf("	mulf %s, %s, %s\n", fregister_str[src_register1], fregister_str[src_register2],
				fregister_str[dest_register]);
			break;
		case 20 :	// divf
			printf("	divf %s, %s, %s\n", fregister_str[src_register1], fregister_str[src_register2],
				fregister_str[dest_register]);
			break;
		case 21 :	// seteqf
			printf("	seteqf %s, %s, %s\n", fregister_str[src_register1], fregister_str[src_register2],
				register_str[dest_register]);
			break;
		case 22 :	// setgrf
			printf("	setgrf %s, %s, %s\n", fregister_str[src_register1], fregister_str[src_register2],
				register_str[dest_register]);
			break;
		case 23 :	// Transcendental functions and square root.
			if (!include_expensive_math) {
				printf("	noop\n");
				break;
			}
			switch ((program[ip] >> 5) & 3) {
			case 0 :
				printf("	sqrtf");
				break;
			case 1:
				printf("	sinf");
				break;
			case 2 :
				printf("	cosf");
				break;
			case 3 :
				printf("	expf");
				break;
			}
			printf(" %s, %s\n", fregister_str[src_register1], fregister_str[dest_register]);
			break;
		case 24 :	// movcf
			printf("	movcf $%8.6ff, %s\n", *(float *)&program[max_instructions + 1 +
				((program[ip] >> 5) & (NU_FLOAT_CONSTANTS - 1))], fregister_str[dest_register]);
			break;
		case 25 :	// cnvrtf
			printf("	cnvrtf %s, %s\n", register_str[src_register1], fregister_str[dest_register]);
			break;
		case 26 :	// roundf
			switch ((program[ip] >> 5) & 3) {
			case 0 :
				printf("	roundzf");
				break;
			case 1 :
				printf("	roundnf");
				break;
			case 2 :
				printf("	roundif");
				break;
			case 3 :
				printf("	roundnif");
				break;
			}
			printf(" %s, %s\n", fregister_str[src_register1], register_str[dest_register]);
			break;
		default :	// noop
			printf("	noop\n");
			break;
		}
		ip++;
		if (ip >= nu_instructions)
			break;
	}
}

static int evaluate_program(unsigned int *program, int *r, float *f) {
	int ip = 0;
	int penalty = 0;	// Keep a penalty score for expensive instructions executed.
	int nu_instructions = (int)program[max_instructions];
	for (;;) {
		int src_register1 = (program[ip] >> 17) & (nu_registers - 1);
		int src_register2 = (program[ip] >> 22) & (nu_registers - 1);
		int dest_register = (program[ip] >> 27) & (nu_registers - 1);
		int opcode = program[ip] & 31;
		if (only_floating_point_math && ((opcode >= 0 && opcode <= 8) || (opcode >= 11 && opcode <= 14)))
			opcode |= 16;
		if (only_integer_math)
			opcode &= 15;
		switch (opcode) {
		case 0 :	// mov
			r[dest_register] = r[src_register1];
			break;
		case 1 :	// add
			r[dest_register] = r[src_register1] + r[src_register2];
			break;
		case 2 :	// sub
			r[dest_register] = r[src_register1] - r[src_register2];
			break;
		case 3 :	// mul
			r[dest_register] = r[src_register1] * r[src_register2];
			penalty++;
			break;
		case 4 :	// div
			if (r[src_register2] == 0 || (r[src_register1] == (- 1 << 31) && r[src_register2] == - 1)) {
				// Ungraceful exit counts as max instructions.
				ip = max_instructions;
				penalty += 100;
				goto abort;
			}
			r[dest_register] = r[src_register1] / r[src_register2];
			penalty += 5;
			break;
		case 5 :	// and
			r[dest_register] = r[src_register1] & r[src_register2];
			break;
		case 6 :	// or
			r[dest_register] = r[src_register1] | r[src_register2];
			break;
		case 7 :	// seteq
			if (r[src_register1] == r[src_register2])
				r[dest_register] = 1;
			else
				r[dest_register] = 0;
			break;
		case 8 :	// setgr
			if (r[src_register1] > r[src_register2])
				r[dest_register] = 1;
			else
				r[dest_register] = 0;
			break;
		case 9 :	// jmpz
			if (r[dest_register] == 0)
				ip += (program[ip] >> 5) & 31;
			break;
		case 10: 	// jmpnz
			if (r[dest_register] != 0)
				ip += (program[ip] >> 5) & 31;
			break;
		case 11 : 	// movsew
			if (program[ip] & (1 << 20))	// Sign bit is set.
				r[dest_register] = (0xFFFF << 16) | ((program[ip] >> 5) & 0xFFFF);
			else
				r[dest_register] = (program[ip] >> 5) & 0xFFFF;
			break;
		case 12 : 	// movlo
			r[dest_register] = (r[dest_register] & 0xFFFF0000) | ((program[ip] >> 5) & 0xFFFF);
			break;
		case 13 : 	// movhi
			r[dest_register] = (r[dest_register] & 0xFFFF) | (((program[ip] >> 5) & 0xFFFF) << 16);
			break;
		case 14 :	// movseb
			if (program[ip] & (1 << 12))	// Sign bit is set.
				r[dest_register] = (0xFFFFFF << 8) | ((program[ip] >> 5) & 0xFF);
			else
				r[dest_register] = (program[ip] >> 5) & 0xFF;
			break;
		case 15 :	// end
			return ip + penalty;
		// Floating point instructions.
		case 16 :	// movf
			f[dest_register] = f[src_register1];
			break;
		case 17 :	// addf
			f[dest_register] = f[src_register1] + f[src_register2];
			break;
		case 18 :	// subf
			f[dest_register] = f[src_register1] - f[src_register2];
			break;
		case 19 :	// mulf
			f[dest_register] = f[src_register1] * f[src_register2];
			break;
		case 20 :	// divf
			if (f[src_register2] == 0) {
				// Ungraceful exit counts as max instructions.
				ip = max_instructions;
				penalty += 100;
				goto abort;
			}
			f[dest_register] = f[src_register1] / f[src_register2];
			penalty += 2;
			break;
		case 21 :	// seteqf
			if (f[src_register1] == f[src_register2])
				r[dest_register] = 1;
			else
				r[dest_register] = 0;
			break;
		case 22 :	// setgrf
			if (f[src_register1] > f[src_register2])
				r[dest_register] = 1;
			else
				r[dest_register] = 0;
			break;
		case 23 :	// Transcendentals and square root.
			if (!include_expensive_math)
				break;	// No-op.
			switch ((program[ip] >> 5) & 3) {
			case 0 :
				if (f[src_register1] < 0) {
					ip = max_instructions;
					penalty += 100;
					goto abort;
				}
				f[dest_register] = sqrtf(f[src_register1]);
				penalty += 10;
				break;
			case 1 :
				f[dest_register] = sinf(f[src_register1]);
				penalty += 10;
				break;
			case 2 :
				f[dest_register] = cosf(f[src_register1]);
				penalty += 10;
				break;
			case 3 :
				f[dest_register] = expf(f[src_register1]);
				penalty += 10;
				break;
			}
			break;
		case 24 :	// movcf
			// Move floating point constant to float register.
			f[dest_register] = *(float *)&program[max_instructions + 1 + ((program[ip] >> 5) &
				(NU_FLOAT_CONSTANTS - 1))];
			break;
		case 25 : 	// cnvrtf
			f[dest_register] = r[src_register1];
			break;
		case 26 :	// roundf
			r[dest_register] = custom_round(f[src_register1], (program[ip] >> 5) & 3);
			break;
		default :	// noop
			break;
		}
		ip++;
		if (ip >= nu_instructions)
			break;
	}
abort :
	if (ip > nu_instructions)
		// Add extra penalty for jumping too far.
		return nu_instructions + penalty + 1;
	else
		return ip + penalty;
}

static double gp_calculate_fitness(const FgenPopulation *pop, const unsigned char *bitstring) {
	int r[32];	// Integer registers
	float f[32];	// Floating-point registers
	double fitness = 0;
	int model_input_correct_count = 0;
	int max_program_cost = INT_MIN;
	for (int i = 0; i < nu_model_inputs; i++) {
		for (int j = 0; j < nu_input_registers; j++)
			r[j] = model_input[i][j];
		// Initialize unused registers with zeroes.
		for (int j = nu_input_registers; j < nu_registers; j++)
			r[j] = 0;
		for (int j = 0; j < nu_registers; j++)
			f[j] = 0;
		int n = evaluate_program((unsigned int *)bitstring, r, f);
		if (n > max_program_cost)
			max_program_cost = n;
		int count = 0;
		for (int j = 0; j < nu_output_registers; j++)
			if (r[j] == model_output[i][j]) {
				fitness += 1.0;
				count++;
			}
		if (count == nu_output_registers)
			model_input_correct_count++;
	}
	if (model_input_correct_count == nu_model_inputs)
		// If the program output was 100% correct, also weigh the cost of the program.
		if (max_instructions - max_program_cost < - 90)
			fitness += 10;
		else
			fitness += 100 + max_instructions - max_program_cost;
	return fitness;
}

// Floating point model version uses least squares technique to arrive at fitness.

static double gp_calculate_fitness_float_model(const FgenPopulation *pop, const unsigned char *bitstring) {
	int r[32];	// Integer registers
	float f[32];	// Floating-point registers
	int max_program_cost = INT_MIN;
	double error = 0;
	for (int i = 0; i < nu_model_inputs; i++) {
		for (int j = 0; j < nu_input_registers; j++)
			f[j] = fmodel_input[i][j];
		for (int j = nu_input_registers; j < nu_registers; j++)
			f[j] = 0;
		for (int j = 0; j < nu_registers; j++)
			r[j] = 0;
		int n = evaluate_program((unsigned int *)bitstring, r, f);
		if (n > max_program_cost)
			max_program_cost = n;
		for (int j = 0; j < nu_output_registers; j++)
			error += (f[j] - fmodel_output[i][j]) * (f[j] - fmodel_output[i][j]);
		// Avoid occurence of NaN in fitness.
		if (isnan(error))
			error = INFINITY;
	}
	double fitness = 1000 / (error + 1);
	if (error < 0.01)
		// If the program output was almost 100% correct, also weigh the cost of the program.
		if (max_instructions - max_program_cost < - 90)
			fitness += 10;
		else
			fitness += 100 + max_instructions - max_program_cost;
	return fitness;
}

static void gp_generation_callback(FgenPopulation *pop, int generation) {
	if (generation % 100 == 0) {
		FgenIndividual *ind = fgen_best_individual_of_population(pop);
		printf("Gen = %d Best fitness = %lf\n", generation, ind->fitness);
		if (fgen_is_cached(pop))
			printf("Fitness cache hit-rate = %lf.\n", fgen_get_cache_hit_rate(pop));
		print_program((unsigned int *)ind->bitstring);
	}
//	if (generation == 10000) {
//		fgen_signal_stop(pop);
//	}
}

// Define the range of the floating point constants used in the chromosome.

#define DOMAIN_MIN (- 100.0)
#define DOMAIN_MAX 100.0

// Custom seed operator because we have a value holding the number of instructions and floating point constants
// at the end of the instruction list.

static void gp_seed(FgenPopulation *pop, unsigned char *bitstring) {
	FgenRNG *rng = fgen_get_rng(pop);
	// Seed with a full chromosome of instructions.
	for (int i = 0; i < max_instructions; i++)
		*(unsigned int *)&bitstring[i * 4] = fgen_random_32(rng);
	*(int *)&bitstring[max_instructions * 4] = max_instructions;
	for (int i = max_instructions + 1; i < max_instructions + 1 + NU_FLOAT_CONSTANTS; i++)
		*(float *)&bitstring[i * 4] = fgen_random_f(rng, DOMAIN_MAX - DOMAIN_MIN) + DOMAIN_MIN;
}

// Custom mutations for the instruction list.

static void gp_mutation_delete(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	unsigned int *p = (unsigned int *)parent;
	unsigned int *c = (unsigned int *)child;
	FgenRNG *rng = fgen_get_rng(pop);
	int nu_instructions = *(int *)&p[max_instructions];
	if (nu_instructions == 0)
		return;
	int i = fgen_random_n(rng, nu_instructions);
	if (i != nu_instructions - 1)
		memmove(&c[i], &c[i + 1], (nu_instructions - 1 - i) * 4);
	(*(int *)&c[max_instructions])--;
}

static void gp_mutation_insert(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	unsigned int *p = (unsigned int *)parent;
	unsigned int *c = (unsigned int *)child;
	FgenRNG *rng = fgen_get_rng(pop);
	int nu_instructions = *(int *)&p[max_instructions];
	int i;
	if (nu_instructions == 0)
		i = 0;
	else
		i = fgen_random_n(rng, nu_instructions);
	if (i < nu_instructions && nu_instructions < max_instructions)
		memmove(&c[i + 1], &c[i], (nu_instructions - i) * 4);
	c[i] = fgen_random_32(rng);
	if (nu_instructions < max_instructions)
		(*(int *)&c[max_instructions])++;
}

// Because we have extra non-bitstring data in the chromosome, we can't use the built-in bitstring
// mutation functions for bitwise mutation of the instructions.

static void gp_mutation_per_bit(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	FgenRNG *rng = fgen_get_rng(pop);
	for (int i = 0; i < max_instructions * 32; i++) {
		int r = fgen_random_16(rng);
		if (r < pop->mutation_probability) {
			/* Mutate the bit. */
			fgen_mutate_bit(child, i);
		}
	}
}

static void gp_mutation_bitstring(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	gp_mutation_per_bit(pop, parent, child);
	FgenRNG *rng = fgen_get_rng(pop);
	if (pop->macro_mutation_probability > 0) {
		int i;
		/* For each data element (instruction). */
		for (i = 0; i < max_instructions; i++) {
			int r = fgen_random_16(rng);
			if (r < pop->macro_mutation_probability)
				/* Completely replace the instruction. */
				*(unsigned int *)&child[i * 4] = fgen_random_32(rng);
		}
	}
}

// The floating point constants at the end of the chromosome need special mutation operators.

/* Custom mutation function that assigns random values within the range to parameters with a mutation */
/* probability. When the function is called child already contains a copy of parent */

static void gp_mutation_float_constants_large(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	float *params_child = (float *)child;
	int i; 
	/* Mutate each parameter with the given probability. */
	for (i = max_instructions + 1; i < max_instructions + 1 + NU_FLOAT_CONSTANTS; i++)
		if (fgen_random_f(pop->rng, 1.0) < pop->macro_mutation_probability_float)
			params_child[i] = fgen_random_f(pop->rng, DOMAIN_MAX - DOMAIN_MIN) + DOMAIN_MIN;
}

/* Custom mutation function that tweaks the parameter by a maximum of 1/20th of the range. */

static void gp_mutation_float_constants_small(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	float *params_parent = (float *)parent;
	float *params_child = (float *)child;
	int i;
	/* Mutate each parameter with the given probability. */
	for (i = max_instructions + 1; i < max_instructions + 1 + NU_FLOAT_CONSTANTS; i++)
		if (fgen_random_f(pop->rng, 1.0) < 0.03) {
			float tweak = fgen_random_f(pop->rng, (DOMAIN_MAX - DOMAIN_MIN) / 10) -
				(DOMAIN_MAX - DOMAIN_MIN) / 20;
			params_child[i] = params_parent[i] + tweak;
			/* Clamp to range. */
			if (params_child[i] < DOMAIN_MIN)
				params_child[i] = DOMAIN_MIN;
			if (params_child[i] > DOMAIN_MAX)
				params_child[i] = DOMAIN_MAX;
		}
}

static void gp_mutation_float_constants(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	gp_mutation_float_constants_small(pop, parent, child);
	gp_mutation_float_constants_large(pop, parent, child);
}

static void gp_mutation(FgenPopulation *pop, const unsigned char *parent, unsigned char *child) {
	FgenRNG *rng = fgen_get_rng(pop);
	float r = fgen_random_f(rng, 1.0);
	if (r < 0.1) {
		gp_mutation_delete(pop, parent, child);
	}
	else
	if (r < 0.2) {
		gp_mutation_insert(pop, parent, child);
	}
	else
		gp_mutation_bitstring(pop, parent, child);
	if (!only_integer_math)
		gp_mutation_float_constants(pop, parent, child);
	return;
}

// Customer crossover function. Exchange up to four contiguous instructions from random locations in each parent.
// The size of the segment is the same for each parent.

static void gp_crossover(FgenPopulation *pop, const unsigned char *parent1, const unsigned char *parent2,
unsigned char *child1, unsigned char *child2) {
	// Copy the whole chromosomes first.
	memcpy(child1, parent1, (max_instructions + 1 + NU_FLOAT_CONSTANTS) * 4);
	memcpy(child2, parent2, (max_instructions + 1 + NU_FLOAT_CONSTANTS) * 4);
	FgenRNG *rng = fgen_get_rng(pop);
	int nu_instructions1 = *(int *)&parent1[max_instructions * 4];
	if (nu_instructions1 == 0)
		return;
	int index1 = fgen_random_n(rng, nu_instructions1);
	int max_size = nu_instructions1 - index1;
	if (max_size > 4)
		max_size = 4;
	int nu_instructions2 = *(int *)&parent2[max_instructions * 4];
	if (nu_instructions2 == 0)
		return;
	if (max_size > nu_instructions2)
		max_size = nu_instructions2;
	int size = fgen_random_n(rng, max_size) + 1;
	int index2 = fgen_random_n(rng, nu_instructions2 - size + 1);
	// Copy the segments.
	memcpy(&child2[index2 * 4], &parent1[index1 * 4], size * 4);
	memcpy(&child1[index1 * 4], &parent2[index2 * 4], size * 4);
}

void run_ga() {
	FgenPopulation *pop;
	FgenCalculateFitnessFunc calculate_fitness_func;
	if (model_floating_point_mode)
		calculate_fitness_func = gp_calculate_fitness_float_model;
	else
		calculate_fitness_func = gp_calculate_fitness;
	// Chromosome:
	// 1. max_instructions x 32-bit instructions.
	// 2. 32-bit number of actual instructions.
	// 3. NU_FLOAT_CONSTANTS floats.
	pop = fgen_create(
		1024,							/* Population size. */
		max_instructions * 32 + 32 + NU_FLOAT_CONSTANTS * 32,	/* Individual size in bits. */
		32,							/* Data element size. */
		gp_generation_callback,
		calculate_fitness_func,
		gp_seed,
		gp_mutation,
		gp_crossover
//		fgen_crossover_uniform_per_element
//		fgen_crossover_noop
		);
	fgen_set_parameters(
		pop,
		FGEN_ELITIST_SUS,		/* Selection type. */
		FGEN_SUBTRACT_MIN_FITNESS,	/* Selection fitness type. */
		0.200,				/* Crossover probability. */
		0.015,				/* Mutation probability per bit. */
		0.05);				/* Macro-mutation probability. */
	fgen_random_seed_with_timer(fgen_get_rng(pop));
	fgen_run_threaded(pop, - 1);
}

