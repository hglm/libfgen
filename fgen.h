/*
    fgen.h -- main API header file.

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


/*
 * User API for the libfgen library.
 */

#ifndef __FGEN_H__
#define __FGEN_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS /* empty */
#define __END_DECLS /* empty */
#endif

// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define FGEN_HELPER_SHARED_IMPORT __declspec(dllimport)
  #define FGEN_HELPER_SHARED_EXPORT __declspec(dllexport)
  #define FGEN_HELPER_SHARED_LOCAL
#else
  #if __GNUC__ >= 4
    #define FGEN_HELPER_SHARED_IMPORT __attribute__ ((visibility ("default")))
    #define FGEN_HELPER_SHARED_EXPORT __attribute__ ((visibility ("default")))
    #define FGEN_HELPER_SHARED_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define FGEN_HELPER_SHARED_IMPORT
    #define FGEN_HELPER_SHARED_EXPORT
    #define FGEN_HELPER_SHARED_LOCAL
  #endif
#endif

// Now we use the generic helper definitions above to define FGEN_API and FGEN_LOCAL.
// FGEN_API is used for the public API symbols. It either DLL imports or DLL exports (or does nothing
// for static build). FGEN_LOCAL is used for non-api symbols.

#ifdef FGEN_SHARED
  // Defined if FGEN is compiled as a shared library.
  #ifdef FGEN_SHARED_EXPORTS
    // Defined if we are building the libfgen shared library (instead of using it).
    #define FGEN_API FGEN_HELPER_SHARED_EXPORT
  #else
    #define FGEN_API FGEN_HELPER_SHARED_IMPORT
  #endif // FGEN_SHARED_EXPORTS
  #define FGEN_LOCAL FGEN_HELPER_SHARED_LOCAL
#else
  // FGEN_SHARED is not defined: this means libfgen is a static lib.
  #define FGEN_API
  #define FGEN_LOCAL
#endif // FGEN_SHARED

#include <stdint.h>  // uint64_t is used.

__BEGIN_DECLS

/*
 * fgen functionality.
 */

/** \defgroup group_fgen fgen (Genetic Algorithm) API
 * @{
 */

typedef struct FGEN_API {
	unsigned char *bitstring;
	double fitness;
	int refcount;
	/* Flags: */
	unsigned char fitness_is_valid;
	unsigned char obsolete;
	unsigned char is_elite;
        unsigned char unused;
} FgenIndividual;

typedef struct FgenPopulation_t FgenPopulation;

typedef struct FGEN_API {
	unsigned char *bitstring;
	double fitness;
	int date_mru;
} FgenCacheEntry;

typedef struct FGEN_API {
	int size;
	int nu_accesses;
	int nu_hits;
	double hit_rate;
	FgenCacheEntry **entry;
	int refcount;
} FgenCache;

#define FGEN_RNG_STATE_SIZE 8
#define FGEN_RNG_STORAGE_SIZE 64

typedef struct FGEN_API {
	unsigned int Q[FGEN_RNG_STATE_SIZE];
	unsigned int c;
	int index;
#if FGEN_RNG_STORAGE_SIZE == 64
	// 64 bits of storage.
	uint64_t storage;
#else
	unsigned int storage;
#endif
	int storage_size;
	unsigned int last_power_of_two;
        unsigned int last_general_range;
	// The bit shift corresponding to the last power of two (log2(n)).
	unsigned char last_power_of_two_shift;
        unsigned char last_general_range_shift;
} FgenRNG;

/** The generation callback function. The current population and the current generation are passed as arguments. */
typedef void (*FgenGenerationCallbackFunc)(FgenPopulation *pop, int generation);
/** The fitness evalution function. The current population and the bitstring to be evaluated are passed as
 * arguments. */
typedef double (*FgenCalculateFitnessFunc)(const FgenPopulation *pop, const unsigned char *bitstring);
/** The seeding operator function type. bitstring is already allocated. */
typedef void (*FgenSeedFunc)(FgenPopulation *pop, unsigned char *bitstring);
/** The mutation operator function type. When this is called the child is allocated and already contains
 * a copy of the parent. This function is called for every individual of the population. */
typedef void (*FgenMutationFunc)(FgenPopulation *pop, const unsigned char *parent,
unsigned char *child);
/** The crossover operator function type. The children are allocated but not filled in. This function is called
  * depending on the crossover probability. */
typedef void (*FgenCrossoverFunc)(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2);

struct FgenPopulation_t {
/* Core parameters. */
	int size;			/**< Size of the population. */
	int individual_size_in_bits;	/**< Size of the bitstring of an individual in bits (multiple of 8). */
	int data_element_size;		/**< Size of a date element within a bitstring (multiple of 8 or 1). */
	int generation;			/**< The current generation. */
	FgenIndividual **ind;		/**< Array of pointers to individuals. */
	int island;			/**< The island to which this population belongs. */
/* Evolution parameters. */
	int selection_type;		/**< The selection type. */
	int tournament_size;		/**< The tournament size for tournament selection. */
	int nu_elites;			/**< The number of elites for elitist selection. */
	int selection_fitness_type;	/**< The fitness selection type. */
	float crossover_probability_float;	/**< The crossover probability. */
	int crossover_probability;
	float mutation_probability_float;	/**< The mutation probability. */
	int mutation_probability;
	float macro_mutation_probability_float;	/**< The macro mutation probability. */
	int macro_mutation_probability;
	float migration_probability_float;	/**< The migration probability per individual. */
	int migration_probability;
	int migration_interval;			/**< The migration interval. */
	int permutation_size;			/**< The permutation size for permutation operators. */
	int generation_callback_interval;	/**< The interval at which the generation callback function is called. */
/* Miscellaneous. */
	int stop_signalled;			/**< Whether a stop is signalled. */
	int cache_is_shared;			/**< Obsolete. */
	int max_threads;			/**< Maximum number of concurrent threads to use. */
	FgenCache *cache;			/**< Pointer to the cache. */
	FgenRNG *rng;				/**< Pointer to random number generator. */
	void *user_data;			/**< User data, used by ffit. */
/* Callback functions. */
	FgenGenerationCallbackFunc fgen_generation_callback_func;	/**< Generation callback function. */
	FgenCalculateFitnessFunc fgen_calculate_fitness_func;		/**< Fitness calculation function. */
	FgenSeedFunc fgen_seed_func;		/**< Seeding operator. */
	FgenMutationFunc fgen_mutation_func;	/**< Mutation operator. */
	FgenCrossoverFunc fgen_crossover_func;	/**< Crossover operator. */
	float *fast_mutation_cumulative_chance;	/**< Used internally. */
	int fast_mutation_nu_bits_to_mutate;	/**< Used internally. */
	int initialization_type;		/**< Flag that can be used to continue with an existing population. */
	float fast_mutation_probability;	/**< Used internally. */
        int nu_data_elements;			/**< Used internally (number of data elements in an individual). */
        int population_size_shift;		/**< Used internally (>= 0 if population size is power of two). */
        int individual_size_shift;		/**< Used internally (>= 0 if individual size is power of two). */
        int data_element_size_shift;		/**< Used internally (>= 0 if data element size is power of two). */
};

/* Fitness selection types. */
#define FGEN_FITNESS_PROPORTIONAL		1
#define FGEN_SUBTRACT_MIN_FITNESS		2
#define FGEN_SUBTRACT_MIN_FITNESS_DIV_2	3

/* Selection types. */
#define FGEN_STOCHASTIC_TYPE_MASK	3
#define FGEN_ELITIST_ELEMENT		4
#define FGEN_EXTINCTION_ELEMENT		8
#define FGEN_STOCHASTIC			0
#define FGEN_SUS			1
#define FGEN_TOURNAMENT			2
#define FGEN_RANK			3
#define FGEN_KILL_TOURNAMENT_ELEMENT	16
#define FGEN_ELITIST_STOCHASTIC			(FGEN_STOCHASTIC | FGEN_ELITIST_ELEMENT)
#define FGEN_ELITIST_SUS			(FGEN_SUS | FGEN_ELITIST_ELEMENT)
#define FGEN_ELITIST_SUS_WITH_EXTINCTION 	(FGEN_SUS | FGEN_ELITIST_ELEMENT | FGEN_EXTINCTION_ELEMENT)		
#define FGEN_ELITIST_TOURNAMENT			(FGEN_TOURNAMENT | FGEN_ELITIST_ELEMENT)
#define FGEN_ELITIST_TOURNAMENT_WITH_EXTINCTION (FGEN_TOURNAMENT | FGEN_ELITIST_ELEMENT | FGEN_EXTINCTION_ELEMENT)
#define FGEN_ELITIST_RANK			(FGEN_RANK | FGEN_ELITIST_ELEMENT)
#define FGEN_KILL_TOURNAMENT			(FGEN_TOURNAMENT | FGEN_KILL_TOURNAMENT_ELEMENT)

/* Initialization types. */
#define FGEN_INITIALIZATION_SEED	0
#define FGEN_INITIALIZATION_CONTINUE	1

/* Flags (used internally). */
#define FGEN_FLAG_POPULATION_SIZE_POWER_OF_TWO	1
#define FGEN_FLAG_INDIVIDUAL_SIZE_POWER_OF_TWO	2

/* Main interface. */

FGEN_API FgenPopulation *fgen_create(
	int population_size,
	int individual_size_in_bits,
	int data_element_size,
	FgenGenerationCallbackFunc fgen_generation_callback_func,
	FgenCalculateFitnessFunc fgen_calculate_fitness_func,
	FgenSeedFunc fgen_seed_func,
	FgenMutationFunc fgen_mutation_func,
	FgenCrossoverFunc fgen_crossover_func
	);
FGEN_API void fgen_initialize(
	FgenPopulation *pop,
	int population_size,
	int individual_size_in_bits,
	int data_element_size,
	FgenGenerationCallbackFunc fgen_generation_callback_func,
	FgenCalculateFitnessFunc fgen_calculate_fitness_func,
	FgenSeedFunc fgen_seed_func,
	FgenMutationFunc fgen_mutation_func,
	FgenCrossoverFunc fgen_crossover_func
	);
FGEN_API void fgen_set_parameters(
	FgenPopulation *pop,
	int selection_type,
	int selection_fitness_type,
	float crossover_probability_float,
	float mutation_per_bit_probability_float,
	float macro_mutation_probability_float
	);
FGEN_API void fgen_run(FgenPopulation *pop, int max_generation);
FGEN_API void fgen_run_threaded(FgenPopulation *pop, int max_generation);
FGEN_API void fgen_destroy(FgenPopulation *pop);
FGEN_API void fgen_run_archipelago(int nu_pops, FgenPopulation **pops, int max_generation);
FGEN_API void fgen_run_archipelago_threaded(int nu_pops, FgenPopulation **pops, int max_generation);
FGEN_API void fgen_run_steady_state(FgenPopulation *pop, int max_generation);
FGEN_API void fgen_run_steady_state_archipelago(int nu_pops, FgenPopulation **pops, int max_generation);
FGEN_API void fgen_run_steady_state_archipelago_threaded(int nu_pops, FgenPopulation **pops, int max_generation);

/* Parameter interface.  */

FGEN_API void fgen_set_mutation_probability(FgenPopulation *pop, float);
FGEN_API void fgen_set_macro_mutation_probability(FgenPopulation *pop, float);
FGEN_API void fgen_set_crossover_probability(FgenPopulation *pop, float);
FGEN_API void fgen_set_selection_fitness_type(FgenPopulation *pop, int);
FGEN_API void fgen_set_selection_type(FgenPopulation *pop, int);
FGEN_API void fgen_set_tournament_size(FgenPopulation *pop, int);
FGEN_API void fgen_set_data_element_size(FgenPopulation *pop, int);
FGEN_API void fgen_set_number_of_elites(FgenPopulation *pop, int);
FGEN_API void fgen_set_permutation_size(FgenPopulation *pop, int);
FGEN_API void fgen_set_user_data(FgenPopulation *pop, void *user_data);
FGEN_API void fgen_enable_cache(FgenPopulation *pop, int cache_size);
FGEN_API void fgen_enable_cache_on_archipelago(int nu_pops, FgenPopulation **pops, int cache_size);
FGEN_API void fgen_invalidate_cache(FgenPopulation *pop);
FGEN_API void fgen_set_migration_probability(FgenPopulation *pop, float);
FGEN_API void fgen_set_migration_interval(FgenPopulation *pop, int);
FGEN_API void fgen_set_generation_callback_interval(FgenPopulation *pop, int);
FGEN_API void fgen_set_initialization_type(FgenPopulation *pop, int);

/* Default operator functions. */

FGEN_API void fgen_seed_random(FgenPopulation *pop, unsigned char *bitstring);
FGEN_API void fgen_mutation_per_bit(FgenPopulation *pop, const unsigned char *parent, unsigned char *child);
FGEN_API void fgen_mutation_per_bit_plus_macro_mutation(FgenPopulation *pop, const unsigned char *parent,
unsigned char *child);
FGEN_API void fgen_mutation_per_bit_fast(FgenPopulation *pop, const unsigned char *parent, unsigned char *child);
FGEN_API void fgen_mutation_per_bit_plus_macro_mutation_fast(FgenPopulation *pop, const unsigned char *parent,
unsigned char *child);
FGEN_API void fgen_crossover_one_point_per_bit(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2);
FGEN_API void fgen_crossover_one_point_per_element(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2);
FGEN_API void fgen_crossover_two_point_per_bit(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2);
FGEN_API void fgen_crossover_two_point_per_element(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2,unsigned char *child1, unsigned char *child2);
FGEN_API void fgen_crossover_uniform_per_bit(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2);
FGEN_API void fgen_crossover_uniform_per_element(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2);
FGEN_API void fgen_seed_permutation_random(FgenPopulation *pop, unsigned char *bitstring);
FGEN_API void fgen_mutation_permutation_swap(FgenPopulation *pop, const unsigned char *parent,
unsigned char *child);
FGEN_API void fgen_mutation_permutation_insert(FgenPopulation *pop, const unsigned char *parent,
unsigned char *child);
FGEN_API void fgen_mutation_permutation_invert(FgenPopulation *pop, const unsigned char *parent,
unsigned char *child);
FGEN_API void fgen_crossover_permutation_order_based(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2);
FGEN_API void fgen_crossover_permutation_position_based(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2);
FGEN_API void fgen_crossover_noop(FgenPopulation *pop, const unsigned char *parent1,
const unsigned char *parent2, unsigned char *child1, unsigned char *child2);

/* Helper functions to be called from the GenerationCallbackFunc. */

FGEN_API FgenIndividual *fgen_best_individual_of_population(FgenPopulation *pop);
FGEN_API FgenIndividual *fgen_worst_individual_of_population(FgenPopulation *pop);
FGEN_API FgenIndividual *fgen_best_individual_of_archipelago(int nu_pops, FgenPopulation **pops);
FGEN_API FgenIndividual *fgen_best_individual_and_island_of_archipelago(int nu_pops, FgenPopulation **pops,
int *island);
FGEN_API void fgen_update_population_fitness(FgenPopulation *pop);
FGEN_API void fgen_invalidate_population_fitness(FgenPopulation *pop);
FGEN_API void fgen_signal_stop(FgenPopulation *pop);
FGEN_API int fgen_individual_size_in_bytes(const FgenPopulation *pop);
FGEN_API int fgen_get_island(const FgenPopulation *pop);
FGEN_API int fgen_is_cached(const FgenPopulation *pop);
FGEN_API double fgen_get_cache_hit_rate(const FgenPopulation *pop);
FGEN_API FgenRNG *fgen_get_rng(const FgenPopulation *pop);
FGEN_API int fgen_get_generation(const FgenPopulation *pop);

/* Bitstring decoding helper functions. */

FGEN_API double fgen_bitstring_uint16_to_double(const unsigned char *bitstring, double domain_min,
double domain_max);
FGEN_API double fgen_bitstring_uint32_to_double(const unsigned char *bitstring, double domain_min,
double domain_max);
FGEN_API double fgen_bitstring_uint64_to_double(const unsigned char *bitstring, double domain_min,
double domain_max);
FGEN_API void fgen_decode_from_gray(const FgenPopulation *pop, const unsigned char *src_bitstring, unsigned char *dest_bitstring);
FGEN_API void fgen_encode_to_gray(const FgenPopulation *pop, const unsigned char *src_bitstring, unsigned char *dest_bitstring);

/* Bitstring helper functions. */

FGEN_API void fgen_mutate_bit(unsigned char *bitstring, int n);
FGEN_API void fgen_set_random_bitstring(FgenRNG *rng, unsigned char *bitstring, int offset, int nu_bits);
FGEN_API int fgen_get_bit(const unsigned char *bitstring, int n);
FGEN_API void fgen_copy_partial_bitstring(const unsigned char *src_bitstring, unsigned char *dest_bitstring, int offset,
int nu_bits);

/* Random number functions. */

FGEN_API FgenRNG *fgen_random_create_rng();
FGEN_API void fgen_random_destroy_rng(FgenRNG *rng);
FGEN_API void fgen_random_seed_rng(FgenRNG *rng, unsigned int seed);
FGEN_API void fgen_random_seed_with_timer(FgenRNG *rng);
FGEN_API int fgen_random_2(FgenRNG *rng);
FGEN_API int fgen_random_8(FgenRNG *rng);
FGEN_API int fgen_random_16(FgenRNG *rng);
FGEN_API unsigned int fgen_random_32(FgenRNG *rng);
FGEN_API unsigned int fgen_random_bits(FgenRNG *rng, int n_bits);
FGEN_API int fgen_random_n(FgenRNG *rng, int n);
// Integer extensions provided with libfgen v0.2.
FGEN_API unsigned int fgen_random_n_max_65536(FgenRNG *rng, unsigned int n);
FGEN_API unsigned int fgen_random_n_max_256(FgenRNG *rng, unsigned int n);
FGEN_API unsigned int fgen_random_n_power_of_two(FgenRNG *rng, unsigned int n);
FGEN_API unsigned int fgen_random_n_power_of_two_max_65536(FgenRNG *rng, unsigned int n);
FGEN_API unsigned int fgen_random_n_power_of_two_max_256(FgenRNG *rng, unsigned int n);
FGEN_API unsigned int fgen_random_n_power_of_two_with_shift(FgenRNG *rng, unsigned int shift);
FGEN_API unsigned int fgen_random_n_power_of_two_repeat(FgenRNG *rng);
FGEN_API unsigned int fgen_random_n_general_repeat(FgenRNG *rng);
FGEN_API void fgen_random_n_power_of_two_prepare_for_repeat(FgenRNG *rng, unsigned int n);
// Prepare for random range that is not necessarily a power two (but allowed to be),
// to be followed by calls to fgen_random_n_general_repeat().
FGEN_API void fgen_random_n_general_prepare_for_repeat(FgenRNG *rng, unsigned int n);
// Efficiently calculate log2(n), returns - 1 if n is not a power of two.
FGEN_API int fgen_calculate_shift(unsigned int n);
// Note: libfgen v0.2 provides significantly higher precision floats and doubles using the
// following standard functions.
FGEN_API float fgen_random_f(FgenRNG *rng, float range);
FGEN_API double fgen_random_d(FgenRNG *rng, double range);
FGEN_API double fgen_random_from_range_d(FgenRNG *rng, double min_bound, double max_bound);
// Floating point extensions provided with libfgen v0.2.
FGEN_API float fgen_random_f_low_precision(FgenRNG *rng, float range);
FGEN_API float fgen_random_d_low_precision(FgenRNG *rng, float range);
FGEN_API float fgen_random_d_high_precision(FgenRNG *rng, float range);

/** @} */

/* 
 * fpso functionality.
 */

/** \defgroup group_fpso fpso (Particle Swarm Optimization) API
  * @{
  */

typedef struct FGEN_API {
	double *position;
	double *velocity;
	double *best_known_position;
	double best_known_error; 
} FpsoIndividual;

typedef struct FpsoPopulation_t FpsoPopulation;
typedef void (*FpsoGenerationCallbackFunc)(FpsoPopulation *pop, int generation);
typedef double (*FpsoCalculateErrorFunc)(const FpsoPopulation *pop, const double *parameters);

struct FpsoPopulation_t {
	int size;
	int nu_params;
	FpsoIndividual *ind;
	double omega;
	double phi1;
	double phi2;
	int topology;
	int bounding_strategy;
	double *lower_bound;
	double *upper_bound;
	double *best_known_position;
	double best_known_error;
	int stop_signalled;
	FgenRNG *rng;
	void *user_data;
        FpsoGenerationCallbackFunc fpso_generation_callback_func;
	FpsoCalculateErrorFunc fpso_calculate_error_func;
};

#define FPSO_DEFAULT_OMEGA 0.7298
#define FPSO_DEFAULT_PHI1 1.49618
#define FPSO_DEFAULT_PHI2 1.49618

#define FPSO_TOPOLOGY_GBEST 0
#define FPSO_TOPOLOGY_LBEST 1

#define FPSO_BOUND_VELOCITY_ELEMENT 1 
#define FPSO_BOUND_POSITION_ELEMENT 2 
#define FPSO_BOUND_NOTHING 0
#define FPSO_BOUND_VELOCITY 1
#define FPSO_BOUND_POSITION 2
#define FPSO_BOUND_POSITION_AND_VELOCITY 3

FGEN_API FpsoPopulation *fpso_create(
	int population_size,
	int nu_parameters,
	FpsoGenerationCallbackFunc fpso_generation_callback_func,
        FpsoCalculateErrorFunc fpso_calculate_error_func
	);
FGEN_API void fpso_set_parameters(
	FpsoPopulation *pop,
	int topology,
	int bounding_strategy,
	double omega,
	double phi1,
	double phi2
	);
FGEN_API void fpso_run(FpsoPopulation *pop, int max_generation);
FGEN_API void fpso_destroy(FpsoPopulation *pop);
FGEN_API void fpso_signal_stop(FpsoPopulation *pop);

/* Parameter interface. */

FGEN_API void fpso_set_parameter_bounds(FpsoPopulation *pso, int parameter_index, double min_bound, double max_bound);
FGEN_API void fpso_set_topology(FpsoPopulation *pso, int type);
FGEN_API void fpso_set_bounding_strategy(FpsoPopulation *pso, int type);
FGEN_API void fpso_set_user_data(FpsoPopulation *pop, void *user_data);

/* Miscelaneous. */

FGEN_API void fpso_bound_position(const FpsoPopulation *pop, const double *source, double *dest);
FGEN_API double *fpso_get_best_known_position(const FpsoPopulation *pop);
FGEN_API double fpso_get_best_known_error(const FpsoPopulation *pop);
FGEN_API FgenRNG *fpso_get_rng(const FpsoPopulation *pop);

/** @} */

/*
 * Ffit functionality.
 */

/** \defgroup group_ffit ffit (Model fitting / function minimization) API
  * @{
  */

typedef struct Ffit_t Ffit;

typedef double (*FfitCalculateErrorFunc)(const Ffit *fit, const double *param);
typedef void (*FfitGenerationCallbackFunc)(Ffit * fit, int generation,
const double *best_param, double best_error);

/** The main data structure for the ffit functionality. */

struct Ffit_t {
	int nu_params;		/**< Number of parameters in the model. */
	int nu_bits_per_param;	/**< Representation in the genetic algorithm (16, 32 or 64). */
	double *range_min;	/**< Minimum value of each parameter. */
	double *range_max;	/**< Maximum value of each parameter. */
	int *mapping;		/**< Mathematical mapping based on the range [0, 1[. */
	int population_size;	/**< Population size of the optimization algorithm. */
	int optimization_type;	/**< Type of optimization algorithm. */
	int stop_signalled;	/**< Flag indicating whether termination of the algorithm is requested. */
	int model_change_signalled;	/**< Flag indicating whether the model function has changed. */
	int nu_islands;			/**< Internal use. */
	double *best_island_params;	/**< Internal use. */
	double best_island_error;	/**< Internal use. */
	int threading_level;	/**< Whether to use threads. */
	int generation_callback_interval;	/**< How often to call the generation callback function. */
	void *population;	/**< Pointer to optimization algorithm population (internal use). */
	FfitGenerationCallbackFunc ffit_generation_callback_func; /**< Generation call-back function. */
	FfitCalculateErrorFunc ffit_calculate_error_func;	  /**< Error calculating function. */
};

#define FFIT_MAPPING_LINEAR		0
#define FFIT_MAPPING_SQUARE		1
#define FFIT_MAPPING_CUBE		2
#define FFIT_MAPPING_LOG		3
#define FFIT_MAPPING_BINOMIAL_1_TO_5	4
#define FFIT_MAPPING_BINOMIAL_1_TO_7	5

#define FFIT_OPTIMIZATION_FGEN_ELEMENT		1
#define FFIT_OPTIMIZATION_ARCHIPELAGO_ELEMENT	2
#define FFIT_OPTIMIZATION_HILL_CLIMB_ELEMENT	4
#define FFIT_OPTIMIZATION_FPSO_ELEMENT		8
#define FFIT_OPTIMIZATION_FGEN_REAL_VALUED_ELEMENT 16
#define FFIT_OPTIMIZATION_FGEN			FFIT_OPTIMIZATION_FGEN_ELEMENT
#define FFIT_OPTIMIZATION_FGEN_ARCHIPELAGO	(FFIT_OPTIMIZATION_FGEN_ELEMENT | FFIT_OPTIMIZATION_ARCHIPELAGO_ELEMENT)
#define FFIT_OPTIMIZATION_FGEN_WITH_HILL_CLIMB	(FFIT_OPTIMIZATION_FGEN_ELEMENT | FFIT_OPTIMIZATION_HILL_CLIMB_ELEMENT)
#define FFIT_OPTIMIZATION_FGEN_ARCHIPELAGO_WITH_HILL_CLIMB	(FFIT_OPTIMIZATION_FGEN_ELEMENT | \
    FFIT_OPTIMIZATION_ARCHIPELAGO_ELEMENT | FFIT_OPTIMIZATION_HILL_CLIMB_ELEMENT)
#define FFIT_OPTIMIZATION_FPSO			FFIT_OPTIMIZATION_FPSO_ELEMENT
#define FFIT_OPTIMIZATION_FGEN_REAL_VALUED	(FFIT_OPTIMIZATION_FGEN_ELEMENT | \
    FFIT_OPTIMIZATION_FGEN_REAL_VALUED_ELEMENT)
#define FFIT_OPTIMIZATION_FGEN_REAL_VALUED_ARCHIPELAGO (FFIT_OPTIMIZATION_FGEN_ELEMENT | \
FFIT_OPTIMIZATION_FGEN_REAL_VALUED_ELEMENT | FFIT_OPTIMIZATION_ARCHIPELAGO_ELEMENT)

#define FFIT_THREADING_DISABLED	0
#define FFIT_THREADING_ENABLED	1

FGEN_API Ffit *ffit_create(
	int nu_params,
	FfitGenerationCallbackFunc ffit_generation_callback_func,
	FfitCalculateErrorFunc ffit_calculate_error_func
	);
FGEN_API void ffit_set_parameter_range_and_mapping(Ffit *fit, int index, double range_min, double range_max, int mapping);

FGEN_API void ffit_run_fgen_with_defaults(Ffit *fit);
FGEN_API void ffit_run_fgen(Ffit *fit, int population_size, int nu_bits_per_param,
int selection_type, FgenCrossoverFunc crossover, float crossover_probability_float,
float mutation_per_bit_probability_float, float macro_mutation_probability_float);
FGEN_API void ffit_run_fgen_archipelago(Ffit *fit, int nu_islands, int population_size,
int nu_bits_per_param, int selection_type, FgenCrossoverFunc crossover,
float crossover_probability_float, float mutation_probability_float,
float macro_mutation_probability_float, float migration_probability_float,
int migration_interval);
FGEN_API void ffit_run_fgen_real_valued(Ffit *fit, int population_size, int selection_type,
float crossover_probability_float, float mutation_per_bit_probability_float, float macro_mutation_probability_float);
FGEN_API void ffit_run_fgen_real_valued_archipelago(Ffit *fit, int nu_islands, int population_size, int selection_type,
float crossover_probability_float, float mutation_probability_float, float macro_mutation_probability_float,
float migration_probability_float, int migration_interval);
/* FGEN_API void ffit_run_fgen_with_hill_climb(Ffit *fit, int population_size, int nu_bits_per_param, int selection_type,
FgenCrossoverFunc crossover, float crossover_probability_float, float mutation_per_bit_probability_float, float
macro_mutation_probability_float); */
FGEN_API void ffit_run_fpso(Ffit *fit, int population_size, int topology, int bounding_strategy,
double omega, double phi1, double phi2);

FGEN_API void ffit_signal_stop(Ffit *fit);
FGEN_API void ffit_signal_model_change(Ffit *fit);
FGEN_API void ffit_destroy(Ffit *fit);
FGEN_API int ffit_get_population_size(const Ffit *fit);
FGEN_API void ffit_get_individual_params(const Ffit *fit, int index, double *params);
FGEN_API void *ffit_get_population(const Ffit *fit);
FGEN_API void ffit_set_threading(Ffit *fit, int threading_level);
FGEN_API void ffit_set_generation_callback_interval(Ffit *fit, int interval);

/** @} */

__END_DECLS

#endif

