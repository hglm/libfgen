/*
    gp.h -- header file for gp, a linear genetic programming example.

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

#define MAX_INPUT_REGISTERS 32
#define MAX_MODEL_INPUTS 256
#define NU_FLOAT_CONSTANTS 8

void run_ga();

void setup_program_square();
void setup_program_count_registers();
void setup_program_greater_or_equal();
void setup_program_even();
void setup_program_abs();
void setup_program_sort();

void setup_program_formula_float();
void setup_program_round_float();
void setup_program_quartic_float();
void setup_program_transcendental_float();
void setup_program_transcendental2_float();

extern int nu_registers;
extern int max_instructions;
extern int nu_input_registers;
extern int nu_output_registers;
extern int nu_model_inputs;
extern int model_floating_point_mode;
extern int only_floating_point_math;
extern int only_integer_math;
extern int include_expensive_math;
extern int *model_input[MAX_MODEL_INPUTS];
extern int *model_output[MAX_MODEL_INPUTS];
extern float *fmodel_input[MAX_MODEL_INPUTS];
extern float *fmodel_output[MAX_MODEL_INPUTS];

