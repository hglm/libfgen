
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "fgen.h"
#include "win32_compat.h"

static double *LocalBestKnownPosition(FpsoPopulation *pop, int index);

/**
 * Create a PSO population.
 * @param population_size The size of the swarm population.
 * @param nu_parameters The number of parameters per individual.
 * @param fpso_generation_callback_func The generation callback function that is called after every generation of the
 * algorithm.
 * @param fpso_calculate_error_func The error evaluation function that calculates the function value or error of an
 * individual.
 * @return The created population.
 */

FpsoPopulation *fpso_create(int population_size, int nu_parameters, FpsoGenerationCallbackFunc
fpso_generation_callback_func, FpsoCalculateErrorFunc fpso_calculate_error_func) {
	int i;
	FpsoPopulation *pop = (FpsoPopulation *)malloc(sizeof(FpsoPopulation));
	pop->size = population_size;
	pop->nu_params = nu_parameters;
	pop->fpso_generation_callback_func = fpso_generation_callback_func;
	pop->fpso_calculate_error_func = fpso_calculate_error_func;
	pop->ind = (FpsoIndividual *)malloc(sizeof(FpsoIndividual) * population_size);
	for (i = 0; i < population_size; i++) {
		pop->ind[i].position = (double *)malloc(sizeof(double) * nu_parameters);
		pop->ind[i].velocity = (double *)malloc(sizeof(double) * nu_parameters);
		pop->ind[i].best_known_position = (double *)malloc(sizeof(double) * nu_parameters);
        }
	pop->lower_bound = (double *)malloc(sizeof(double) * nu_parameters);
	pop->upper_bound = (double *)malloc(sizeof(double) * nu_parameters);
	pop->best_known_position = (double *)malloc(sizeof(double) * nu_parameters);
	pop->stop_signalled = 0;
	pop->rng = fgen_random_create_rng();
	return pop;
}

/**
 * Set the parameters for an fpso population.
 * @param pop The population.
 * @param topology The swarm topology. One of FPSO_TOPOLOGY_GBEST and FPSO_TOPOLOGY_LBEST.
 * @param bounding_strategy The bounding strategy for parameters used in the algorithm. One of the following:
 * - FPSO_BOUND_NOTHING
 * - FPSO_BOUND_POSITION
 * - FPSO_BOUND_VELOCITY
 * - FPSO_BOUND_POSITION_AND_VELOCITY
 * Irrespective of the value, the FpsoCalculateErrorFunc will always received bounded parameters.
 * @param omega The omega value used in the PSO. A default value is defined as FPSO_DEFAULT_OMEGA.
 * @param phi1 The phi1 value used in the PSO. A default value is defined as FPSO_DEFAULT_PHI1.
 * @param phi2 the phi2 value used in the PSO. A default value is defined as FPSO_DEFAULT_PHI2.
 */

void fpso_set_parameters(FpsoPopulation *pop, int topology, int bounding_strategy, double omega, double phi1,
double phi2) {
	pop->topology = topology;
	pop->bounding_strategy = bounding_strategy;
	pop->omega = omega;
	pop->phi1 = phi1;
	pop->phi2 = phi2;
}

/**
 * Destroy all data structures associated with a PSO population.
 * @param pop The PSO population.
 */

void fpso_destroy(FpsoPopulation *pop) {
	int i;
	fgen_random_destroy_rng(pop->rng);
	for (i = 0; i < pop->size; i++) {
		free(pop->ind[i].position);
		free(pop->ind[i].velocity);
		free(pop->ind[i].best_known_position);
	}
	free(pop->ind);
	free(pop->lower_bound);
	free(pop->upper_bound);
	free(pop->best_known_position);
	free(pop);
}

static double bound(double x, double l, double u) {
	if (x < l)
		return l;
	if (x > u)
		return u;
	return x;
}

/**
 * Bound the parameters given in src to the predefined parameter ranges and store them in dest.
 * @param pop The PSO population.
 * @param src The source parameters.
 * @param dest The location where the bounded parameters are stored.
 */

void fpso_bound_position(const FpsoPopulation *pop, const double *src, double *dest) {
	int i;
	for (i = 0; i < pop->nu_params; i++)
		dest[i] =
		    bound(src[i], pop->lower_bound[i], pop->upper_bound[i]);
}

static void CopyPosition(FpsoPopulation *pop, const double *src, double *dest) {
	memcpy(dest, src, sizeof(double) * pop->nu_params);
}

static void InitializePopulation(FpsoPopulation *pop) {
	int i;
	pop->best_known_error = POSITIVE_INFINITY_DOUBLE;
	for (i = 0; i < pop->size; i++) {
		int j;
		// Initialize position.
		for (j = 0; j < pop->nu_params; j++)
			pop->ind[i].position[j] = fgen_random_from_range_d(pop->rng, pop->lower_bound[j],
				pop->upper_bound[j]);
		CopyPosition(pop, pop->ind[i].position, pop->ind[i].best_known_position);
		pop->ind[i].best_known_error = pop->fpso_calculate_error_func(pop, pop->ind[i].best_known_position);
		if (pop->topology == FPSO_TOPOLOGY_GBEST) {
			// Update the best known position of the population.
			double f = pop->ind[i].best_known_error;
			if (f < pop->best_known_error) {
				pop->best_known_error = f;
				CopyPosition(pop, pop->ind[i].position, pop->best_known_position);
			}
		}
		/* Initialize velocity. */
		for (j = 0; j < pop->nu_params; j++)
			pop->ind[i].velocity[j] = fgen_random_from_range_d(pop->rng,
				- fabs(pop->upper_bound[j] - pop->lower_bound[j]),
				fabs(pop->upper_bound[j] - pop->lower_bound[j]));
	}
}

static void UpdateVelocity(FpsoPopulation *pop, int i, const double *local_best_known_position) {
	int j;
	for (j = 0; j < pop->nu_params; j++) {
		double r_1 = fgen_random_from_range_d(pop->rng, 0, pop->phi1);
		double r_2 = fgen_random_from_range_d(pop->rng, 0, pop->phi2);
		if (pop->topology == FPSO_TOPOLOGY_GBEST)
			pop->ind[i].velocity[j] = pop->omega * pop->ind[i].velocity[j] +
				r_1 * (pop->ind[i].best_known_position[j] - pop->ind[i].position[j]) +
				r_2 * (pop->best_known_position[j] - pop->ind[i].position[j]);
		else	// topology == FPSO_TOPOLOGY_LBEST
			pop->ind[i].velocity[j] = pop->omega * pop->ind[i].velocity[j] +
				r_1 * (pop->ind[i].best_known_position[j] - pop->ind[i].position[j]) +
				r_2 * (local_best_known_position[j] - pop->ind[i].position[j]);
		if (pop->bounding_strategy & FPSO_BOUND_VELOCITY_ELEMENT)
			pop->ind[i].velocity[j] = bound(pop->ind[i].velocity[j],
				-fabs(pop->upper_bound[j] - pop->lower_bound[j]),
				fabs(pop->upper_bound[j] - pop->lower_bound[j]));
	}
}

static void UpdatePosition(FpsoPopulation *pop, int i) {
	int j;
	for (j = 0; j < pop->nu_params; j++) {
		pop->ind[i].position[j] += pop->ind[i].velocity[j];
		if (pop->bounding_strategy & FPSO_BOUND_POSITION_ELEMENT)
			pop->ind[i].position[j] = bound(pop->ind[i].position[j],
				  pop->lower_bound[j], pop->upper_bound[j]);
	}
}

static void UpdateLBESTBestKnownPosition(FpsoPopulation *pop, int i, double *scratch_position) {
	double error_position;
	// If the position is not bound, do bound it for the error calculation only.
	if (!(pop->bounding_strategy & FPSO_BOUND_POSITION_ELEMENT)) {
		fpso_bound_position(pop,pop->ind[i].position, scratch_position);
		error_position = pop->fpso_calculate_error_func(pop, scratch_position);
		if (error_position < pop->ind[i].best_known_error) {
			CopyPosition(pop, scratch_position, pop->ind[i].best_known_position);
			pop->ind[i].best_known_error = error_position;
		}
	} else {
		error_position = pop->fpso_calculate_error_func(pop, pop->ind[i].position);
		if (error_position < pop->ind[i].best_known_error) {
			CopyPosition(pop, pop->ind[i].position, pop->ind[i].best_known_position);
			pop->ind[i].best_known_error = error_position;
		}
	}
}

static void UpdateGBESTBestKnownPositions(FpsoPopulation *pop, double *scratch_position) {
	int i;
	for (i = 0; i < pop->size; i++) {
		double error_position;
		// If the position is not bound, do bound it for the error calculation only.
		if (!(pop->bounding_strategy & FPSO_BOUND_POSITION_ELEMENT)) {
			fpso_bound_position(pop, pop->ind[i].position, scratch_position);
			error_position = pop->fpso_calculate_error_func(pop, scratch_position);
			if (error_position < pop->ind[i].best_known_error) {
				CopyPosition(pop, scratch_position, pop->ind[i].best_known_position);
				pop->ind[i].best_known_error = error_position;
			}
			if (error_position < pop->best_known_error) {
				pop->best_known_error = error_position;
				CopyPosition(pop, scratch_position, pop->best_known_position);
			}
		} else {
			error_position = pop->fpso_calculate_error_func(pop, pop->ind[i].position);
			if (error_position < pop->ind[i].best_known_error) {
				CopyPosition(pop, pop->ind[i].position, pop->ind[i].best_known_position);
				pop->ind[i].best_known_error = error_position;
			}
			if (error_position < pop->best_known_error) {
				pop->best_known_error = error_position;
				CopyPosition(pop, pop->ind[i].position, pop->best_known_position);
			}
		}
	}
}

static void UpdateLBESTGlobalBestKnownPosition(FpsoPopulation *pop, double *scratch_position) {
	int i;
	pop->best_known_error = POSITIVE_INFINITY_DOUBLE;
	for (i = 0; i < pop->size; i++)
		if (pop->ind[i].best_known_error < pop->best_known_error) {
			if (!(pop->bounding_strategy & FPSO_BOUND_POSITION_ELEMENT)) {
//				fpso_bound_position(pop, pop->ind[i].position, scratch_position); Wrong?
				fpso_bound_position(pop, pop->ind[i].best_known_position, scratch_position);
				CopyPosition(pop, scratch_position, pop->best_known_position);
			} else
				CopyPosition(pop, pop->ind[i].best_known_position, pop->best_known_position);
			pop->best_known_error = pop->ind[i].best_known_error;
		}
}

/* Main algorithm. */

/**
 * Run the particle swarm optimization algorithm.
 * @param pop The PSO population.
 * @param max_generation The number of generations to run. The algorithm terminates when this generation is reached. If it
 * is equal to - 1, the algorithm runs indefinitely until fpso_signal_stop() is called.
 */

void fpso_run(FpsoPopulation *pop, int max_generation) {
	double *scratch_position;
	/* Initialize the population. */
	InitializePopulation(pop);

	int generation = 0;
	pop->fpso_generation_callback_func(pop, 0);
	if (pop->stop_signalled)
		return;

	scratch_position = (double *)malloc(sizeof(double) * pop->nu_params);

	/* Run the algorithm. */
	for (;;) {
		int i;
		for (i = 0; i < pop->size; i++) {
			double *local_best_known_position;
			if (pop->topology == FPSO_TOPOLOGY_LBEST)
				local_best_known_position = LocalBestKnownPosition(pop, i);

			UpdateVelocity(pop, i, local_best_known_position);

			UpdatePosition(pop, i);

			/* For LBEST topology, immediately update the best known position. */
			if (pop->topology == FPSO_TOPOLOGY_LBEST)
				UpdateLBESTBestKnownPosition(pop, i, scratch_position);

		}
		/* For GBEST topology, do sychronous updates (update best postions after all positions and */
		/* velocities have been updated. */
		if (pop->topology == FPSO_TOPOLOGY_GBEST)
			UpdateGBESTBestKnownPositions(pop, scratch_position);

		/* For LBEST topology, the best known position of the population is undefined at this stage. */
		/* Before calling the generation_callback_func, update the global best known position. */
		if (pop->topology == FPSO_TOPOLOGY_LBEST)
			UpdateLBESTGlobalBestKnownPosition(pop, scratch_position);

		generation++;
		pop->fpso_generation_callback_func(pop, generation);
		if ((max_generation != - 1 && generation == max_generation) || pop->stop_signalled)
			break;
	}
	free(scratch_position);
}

double *LocalBestKnownPosition(FpsoPopulation *pop, int index) {
/* Ring topology. */
	double *best_position = pop->ind[index].best_known_position;
	double best_error = pop->ind[index].best_known_error;
	/* Check the left neighbour. */
	int index2;
	if (index == 0)
		index2 = pop->size - 1;
	else
		index2 = index - 1;
	if (pop->ind[index2].best_known_error < best_error) {
		best_position = pop->ind[index2].best_known_position;
		best_error = pop->ind[index2].best_known_error;
	}
	/* Check the right neighbour. */
	int index3 = (index + 1) % pop->size;
	if (pop->ind[index3].best_known_error < best_error) {
		best_position = pop->ind[index3].best_known_position;
/*        best_error = ind[index3].best_known_error; */
	}
	return best_position;
}

/* Parameter setting. */

/**
 * Set the bounds for a parameter.
 * @param pop The PSO population.
 * @param index The index (starting with 0) of the parameter of which the bounds are to be set.
 * @param min_bound The lower bound of the range.
 * @param max_bound The upper bound of the range (exclusive).
 */

void fpso_set_parameter_bounds(FpsoPopulation *pop, int index, double min_bound, double max_bound) {
	pop->lower_bound[index] = min_bound;
	pop->upper_bound[index] = max_bound;
}

/**
 * Set the swarm topology.
 * @param pop The PSO population.
 * @param type The swarm topology. One of FPSO_TOPOLOGY_GBEST and FPSO_TOPOLOGY_LBEST.
 */

void fpso_set_topology(FpsoPopulation *pop, int type) {
	pop->topology = type;
}

/**
 * Set the bounding strategy used in the algorithm.
 * @param pop The PSO population.
 * @param type The bounding strategy for parameters used in the algorithm. 
 * Irrespective of the value, the FpsoCalculateErrorFunc will always received bounded parameters.
 * The value must be one of the following:
 * - FPSO_BOUND_NOTHING
 * - FPSO_BOUND_POSITION
 * - FPSO_BOUND_VELOCITY
 * - FPSO_BOUND_POSITION_AND_VELOCITY
 *
 */

void fpso_set_bounding_strategy(FpsoPopulation *pop, int type) {
	pop->bounding_strategy = type;
}


/** 
 * Set the user_data field of the PSO population.
 */

void fpso_set_user_data(FpsoPopulation *pop, void *user_data) {
	pop->user_data = user_data;
}

/**
 * Return the global best known position (parameters with the smallest known error) of the PSO population.
 * @param pop The PSO population.
 * @return A pointer to the best known parameters.
 */

double *fpso_get_best_known_position(const FpsoPopulation *pop) {
	return pop->best_known_position;
}

/**
 * Return the error value associated with the global best known position.
 * @param pop The PSO population.
 * @return The best error value.
 */

double fpso_get_best_known_error(const FpsoPopulation *pop) {
	return pop->best_known_error;
}

/**
 * Signal stop to the PSO algorithm. When called from the generation callback function, the algorithm will
 * terminate after the callback function returns.
 * @param pop The PSO population.
 */

void fpso_signal_stop(FpsoPopulation *pop) {
	pop->stop_signalled = 1;
}

/**
 * Return the random number generator belonging to the PSO population.
 */

FgenRNG *fpso_get_rng(const FpsoPopulation *pop) {
	return pop->rng;
}

