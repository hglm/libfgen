/*
    programs.c -- model programs for gp.

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
#ifdef __GNUC__
#include <fenv.h>
#endif
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

int nu_registers = 8;
int max_instructions = 16;
int nu_input_registers = 1;
int nu_output_registers = 1;
int nu_model_inputs = 1;
int model_floating_point_mode = 0;
int only_floating_point_math = 0;
int only_integer_math = 0;
int include_expensive_math = 0;
int *model_input[MAX_MODEL_INPUTS];
int *model_output[MAX_MODEL_INPUTS];
float *fmodel_input[MAX_MODEL_INPUTS];
float *fmodel_output[MAX_MODEL_INPUTS];

static void allocate_model() {
	if (model_floating_point_mode) {
		for (int i = 0; i < nu_model_inputs; i++) {
			fmodel_input[i] = (float *)malloc(nu_input_registers * sizeof(float));
			fmodel_output[i] = (float *)malloc(nu_output_registers * sizeof(float));
		}
		return;
	}
	for (int i = 0; i < nu_model_inputs; i++) {
		model_input[i] = (int *)malloc(nu_input_registers * sizeof(int));
		model_output[i] = (int *)malloc(nu_output_registers * sizeof(int));
	}
}

// The following program gets passed a random value in r0 should return the value
// r0 * r0 in r1 with r0 set to 0.

void setup_program_square() {
	nu_input_registers = 1;
	nu_output_registers = 2;
	nu_model_inputs = 8;
	max_instructions = 8;
	model_floating_point_mode = 0;
	allocate_model();
	for (int i = 0; i < nu_model_inputs; i++) {
		for (int j = 0; j < nu_input_registers; j++) {
			model_input[i][j] = random() & 255;
		}
		model_output[i][0] = 0;
		model_output[i][1] = model_input[i][0] * model_input[i][0];
	}
}

// The following program has no input values and should return the index number of the register in each of the
// registers.

void setup_program_count_registers() {
	nu_input_registers = 0;
	nu_output_registers = nu_registers;
	nu_model_inputs = 1;
	max_instructions = nu_registers;
	model_floating_point_mode = 0;
	allocate_model();
	for (int i = 0; i < nu_model_inputs; i++) {
		for (int j = 0; j < nu_output_registers; j++) {
			model_output[i][j] = j;
		}
	}
}

// The following program gets passed two random values in r0 and r1 should return 2 in r0 if r0 >= r1 and
// 0 otherwise.

void setup_program_greater_or_equal() {
	nu_input_registers = 2;
	nu_output_registers = 1;
	nu_model_inputs = 64;
	max_instructions = 16;
	model_floating_point_mode = 0;
	only_integer_math = 1;
	allocate_model();
	int f = 0;
	int ok = 0;
again :
	for (int i = 0; i < nu_model_inputs; i++) {
		for (int j = 0; j < nu_input_registers; j++) {
			model_input[i][j] = random() & 31;
		}
		if (model_input[i][0] >= model_input[i][1])
			model_output[i][0] = 2;
		else
			model_output[i][0] = 0;
		// We want at least one model input with equal values.
		if (model_input[i][0] == model_input[i][1])
			ok = 1;
	}
	if (!ok)
		goto again;
}

// The following program gets passed a random value and should return 1 if it is even, 0 otherwise.

void setup_program_even() {
	nu_input_registers = 1;
	nu_output_registers = 1;
	nu_model_inputs = 16;
	max_instructions = 8;
	model_floating_point_mode = 0;
	only_integer_math = 1;
	allocate_model();
	int f = 0;
	int ok = 0;
	for (int i = 0; i < nu_model_inputs; i++) {
		model_input[i][0] = random() & 31;
		if (model_input[i][0] & 1)
			model_output[i][0] = 0;
		else
			model_output[i][0] = 1;
	}
}

// The following program gets passed a random value and should return the absolute value.

void setup_program_abs() {
	nu_input_registers = 1;
	nu_output_registers = 1;
	nu_model_inputs = 16;
	max_instructions = 8;
	model_floating_point_mode = 0;
	only_integer_math = 1;
	allocate_model();
	int f = 0;
	int ok = 0;
	for (int i = 0; i < nu_model_inputs; i++) {
		model_input[i][0] = (random() & 31) - 15;
		if (model_input[i][0] >= 0)
			model_output[i][0] = model_input[i][0];
		else
			model_output[i][0] = - model_input[i][0];
	}
}

// The following program gets passed three integer values in r0, r1 and r3 and should return them
// sorted in increasing order in r0, r1, and r2.

void setup_program_sort() {
	nu_input_registers = 3;
	nu_output_registers = 3;
	nu_model_inputs = 64;
	max_instructions = 24;
	model_floating_point_mode = 0;
	only_integer_math = 1;
	allocate_model();
	for (int i = 0; i < nu_model_inputs; i++) {
		for (int j = 0; j < nu_input_registers; j++) {
			model_input[i][j] = (random() & 0xFF) + 1;
			model_output[i][j] = model_input[i][j];
		}
		if (model_output[i][0] > model_output[i][1]) {
			float temp = model_output[i][0];
			model_output[i][0] = model_output[i][1];
			model_output[i][1] = temp;
		}
		if (model_output[i][1] > model_output[i][2]) {
			float temp = model_output[i][1];
			model_output[i][1] = model_output[i][2];
			model_output[i][2] = temp;
		}
		if (model_output[i][0] > model_output[i][1]) {
			float temp = model_output[i][0];
			model_output[i][0] = model_output[i][1];
			model_output[i][1] = temp;
		}
	}
}

// The following program gets passed a random floating point value in f0 and should return the value
// f0 * f0  + 2.5 * f0 + 5.5 in f1 with f0 unchanged.

void setup_program_formula_float() {
	nu_input_registers = 1;
	nu_output_registers = 2;
	nu_model_inputs = 16;
	max_instructions = 8;
	model_floating_point_mode = 1;
	only_floating_point_math = 1;
	allocate_model();
	for (int i = 0; i < nu_model_inputs; i++) {
		fmodel_input[i][0] = (float)random() * 100.0 / RAND_MAX - 50.0;
		fmodel_output[i][0] = fmodel_input[i][0];
		fmodel_output[i][1] = fmodel_input[i][0] * fmodel_input[i][0] + 2.5 * fmodel_input[i][0] + 5.5;
	}
}

// The following program gets passed a random floating point value in f0 and should return the rounded down
// (integer) value in f0.

void setup_program_round_float() {
	nu_input_registers = 1;
	nu_output_registers = 1;
	nu_model_inputs = 16;
	max_instructions = 8;
	model_floating_point_mode = 1;
	only_floating_point_math = 1;
	allocate_model();
	for (int i = 0; i < nu_model_inputs; i++) {
		fmodel_input[i][0] = (float)random() * 100.0 / RAND_MAX;
		fmodel_output[i][0] = custom_round(fmodel_input[i][0], 3);
	}
}

// The following program gets passed a random floating point value in f0 and should return the value
// f0^4 + f0^3 + f0^2 + f0 in f0.

void setup_program_quartic_float() {
	nu_input_registers = 1;
	nu_output_registers = 1;
	nu_model_inputs = 20;
	max_instructions = 16;
	model_floating_point_mode = 1;
	only_floating_point_math = 1;
	allocate_model();
	for (int i = 0; i < nu_model_inputs; i++) {
		fmodel_input[i][0] = (float)random() * 10.0 / RAND_MAX;
		fmodel_output[i][0] = pow(fmodel_input[i][0], 4) + pow(fmodel_input[i][0], 3) + pow(fmodel_input[i][0], 2)
			+ fmodel_input[i][0];

	}
}

// The following program gets passed a random floating point value in f0 and should return the value
// sin(f0^4 + f0^2) in f0.

void setup_program_transcendental_float() {
	nu_input_registers = 1;
	nu_output_registers = 1;
	nu_model_inputs = 20;
	max_instructions = 10;
	model_floating_point_mode = 1;
	only_floating_point_math = 1;
	include_expensive_math = 1;
	allocate_model();
	for (int i = 0; i < nu_model_inputs; i++) {
		fmodel_input[i][0] = (float)random() * 10.0 / RAND_MAX;
		fmodel_output[i][0] = sinf(pow(fmodel_input[i][0], 4) + pow(fmodel_input[i][0], 2));
	}
}

// The following program gets passed a random floating point value in f0 and should return the value
// sin(exp(sin(exp(sin(f0))))) in f0.

void setup_program_transcendental2_float() {
	nu_input_registers = 1;
	nu_output_registers = 1;
	nu_model_inputs = 20;
	max_instructions = 10;
	model_floating_point_mode = 1;
	only_floating_point_math = 1;
	include_expensive_math = 1;
	allocate_model();
	for (int i = 0; i < nu_model_inputs; i++) {
		fmodel_input[i][0] = (float)random() * 10.0 / RAND_MAX;
		fmodel_output[i][0] = sinf(expf(sinf(expf(sinf(fmodel_input[i][0])))));
	}
}

