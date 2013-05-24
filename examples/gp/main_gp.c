/*
    main.c -- main module of gp, a linear genetic programming example.

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
#include "gp.h"

static int nu_programs = 11;

char *program_name[] = {
	"square", "count_registers", "greater_or_equal", "even", "abs", "sort",
	"formula", "round", "quartic", "transcendental", "transcendental2" };

char *program_description[] = {
	"This program gets passed a random integer value in r0 and returns the squared value in r1 with r0 set to zero.\n",
	"This program has no input values and returns the index number of each register in each integer register.\n",
	"This program gets passed two random integer values in r0 and r1 and returns 2 in r0 if r0 >= r1 and "
	"0 otherwise.\n",
	"This program gets passed a random integer value in r0 and returns 1 in r0 if it is even, 0 otherwise.\n",
	"This program gets passed a random integer value in r0 and returns its absolute value in r0.\n",
	"This program gets passed three random integers in r0, r1 and r2 and sorts them so that r0, r1, r2 "
	"are in increasing order of value.\n",
	"This program gets passed a random floating point value in f0 and returns the result of the formula "
	"f0^2 + 2.5 * f0 + 5.5 in f0.\n",
	"This program gets passed a random floating point value in f0 and returns the rounded down (integer) value in f0.\n",
	"This program gets passed a random floating point value in f0 and returns the result of the formula "
	" f0^4 + f0^3 + f0^2 + f0 in f0.\n",
	"This program gets passed a random floating point value in f0 and returns the value sin(f0^4 + f0^2) in f0.\n",
	"This program gets passed a random floating point value in f0 and returns the value of the formula "
	"sin(exp(sin(exp(sin(f0))))) in f0.\n"
	};

typedef void (*SetupProgramFuncType)();

SetupProgramFuncType setup_program_func[] = {
	setup_program_square,
	setup_program_count_registers,
 	setup_program_greater_or_equal,
	setup_program_even,
	setup_program_abs,
	setup_program_sort,
	// Floating point model programs.
	setup_program_formula_float,
	setup_program_round_float,
	setup_program_quartic_float,
	setup_program_transcendental_float,
	setup_program_transcendental2_float
};

int main(int argc, char **argv) {
	if (argc == 1) {
		printf("gp -- linear genetic programming example using fgen.\n"
			"Usage: gp <options> <program name>\n"
			"Options:\n"
			"-r <number>   Number of registers to use in generated program code (default 8).\n"
			"              Valid values are powers of two up to 32, suggested values are 4,\n"
			"              8, 16 and 32.\n"
			"-d            Display a description of each of the generatable programs.\n"
			"<program name> is the name of the program to be generated, based on a certain\n"
			"number of model inputs and outputs. It must be one from the following list:\n"
			);
		for (int i = 0; i < nu_programs; i++) {
			printf("%s", program_name[i]);
			if (i < nu_programs - 1)
				printf(", ");
		}
		printf("\n");
		exit(0);
	}

	int i = 1;
	for (; i < argc;) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			case 'r' :
				if (i + 1 >= argc) {
					printf("Error -- missing value after option -r.\n");
					exit(1);
				}
				nu_registers = atoi(argv[i + 1]);
				switch (nu_registers) {
				case 1 :
				case 2 :
				case 4 :
				case 8 :
				case 16 :
				case 32 :
					break;
				default :
					printf("Error -- invalid number of registers specified.\n");
					exit(1);
				}
				i += 2;
				break;
			case 'd' :
				for (int j = 0; j < nu_programs; j++)
					printf("%s:\n%s", program_name[j], program_description[j]);
				exit(0);
			}
			continue;
		}
		break;
	}

	if (i >= argc) {
		printf("Must specify program name.\n");
		exit(1);
	}

	int program = - 1;
	for (int j = 0; j < nu_programs; j++)
		if (strcmp(argv[i], program_name[j]) == 0)
			program = j;
	if (program == - 1) {
		printf("Unknown program name. Run gp without arguments for a list of valid names.\n");
		exit(0);
	}

	setup_program_func[program]();
	run_ga();
}

