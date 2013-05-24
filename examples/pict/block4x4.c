/*
    block4x4.c -- texture compression interface of pict, a graphical example using fgen.

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
#include <gtk/gtk.h>
#include <fgen.h>

#include "pict.h"
#include "decode.h"

// Block4x4 texture compression.
//
// External variables:
//
// unsigned char *image_pixels; Pointer to source image pixel data
// int image_width;		The width of the source image.
// int image_height;		The height of the source image.
// int image_extended_width;	The extended width of the source image (width rounded up to multiple of 4).
// int image_extended_height;	The extended height of the source image (height rounded up to a multiple of 4).
// int image_rowstride;		The offset between scanlines of the source image in bytes.
// int nu_generations_per_subdivision;
// int verbose;
//
// copy_buffer_to_result_pixbuf(int sub_image_width, int sub_image_height, int sub_x_offset,
// int sub_y_offset)

typedef struct {
	int etc_modes_allowed;
} PopUserData;

static int sub_x_offset, sub_y_offset;
static unsigned int *sub_image_buffer;

static double block4x4_calculate_fitness(const FgenPopulation *pop, const unsigned char *bitstring) {
	unsigned int image_buffer[16];
	int r;
	int modes_allowed;
	if (pop == NULL)
		modes_allowed = ETC2_MODE_ALLOWED_ALL;
	else
		modes_allowed = ((PopUserData *)pop->user_data)->etc_modes_allowed;
	if (texture_format == 0)
		r = draw_block4x4_etc1(bitstring, image_buffer, modes_allowed);
	else
	if (texture_format == 1)
		r = draw_block4x4_dxt1(bitstring, image_buffer);
	else
		r = draw_block4x4_rgb8_etc2(bitstring, image_buffer, modes_allowed);
	if (r == 0) {
//		printf("Invalid block.\n");
		return 0;	// Fitness is zero for invalid blocks.
	}
	unsigned int *pix = image_buffer;
	unsigned char *pixbuf_pix = image_pixels + sub_y_offset * image_rowstride + sub_x_offset * 3;
	int error = 0;
	for (int y = 0; y < 4; y++) {
		for (int x = 0; x < 4; x++) {
			unsigned char *pixb = (unsigned char *)&pix[x];
			unsigned int r = pixb[0];
			unsigned int g = pixb[1];
			unsigned int b = pixb[2];
			unsigned int r_source = pixbuf_pix[x * 3];
			unsigned int g_source = pixbuf_pix[x * 3 + 1];
			unsigned int b_source = pixbuf_pix[x * 3 + 2];
			error += (r - r_source) * (r - r_source);
			error += (g - g_source) * (g - g_source);
			error += (b - b_source) * (b - b_source);
		}
		pix = &pix[4];
		pixbuf_pix = &pixbuf_pix[image_rowstride];
	}
	return (double)1 / error;
}

static int count = 0;
static unsigned int *compressed_image;

static void block4x4_generation_callback(FgenPopulation *pop, int generation) {
	// Do nothing.
}

static void report_solution(FgenIndividual *best) {
	if (texture_format == 0)
		draw_block4x4_etc1(best->bitstring, sub_image_buffer, ETC_MODE_ALLOWED_ALL);
	else
	if (texture_format == 1)
		draw_block4x4_dxt1(best->bitstring, sub_image_buffer);
	else
		draw_block4x4_rgb8_etc2(best->bitstring, sub_image_buffer, ETC2_MODE_ALLOWED_ALL);
	copy_buffer_to_result_pixbuf(sub_image_buffer, 4, 4, sub_x_offset, sub_y_offset);
//	printf("RMSE of block per pixel: %lf\n", sqrt((1.0 / block4x4_calculate_fitness(NULL, best->bitstring)) / 16));
	// Calculate the block index of the block.
	int compressed_image_index = (sub_y_offset / 4) * (image_extended_width / 4) + sub_x_offset / 4;
	// Copy the 64-bit block.
	compressed_image[compressed_image_index * 2] = *(unsigned int *)best->bitstring;
	compressed_image[compressed_image_index * 2 + 1] = *(unsigned int *)&best->bitstring[4];
	count++;
	if (count == 8 || nu_generations_per_subdivision >= 1000) {
		update_area2();
		count = 0;
	}
	while (gtk_events_pending())
		gtk_main_iteration();
}

static void write_pkm_file() {
	printf("Writing .pkm file.\n");
	FILE *f = fopen(output_filename, "wb");
	if (f == NULL) {
		printf("Error opening output file.\n");
		exit(1);
	}
	fputc('P', f); fputc('K', f); fputc('M', f); fputc(' ', f);
	fputc('1', f); fputc('0', f);
	fputc(0, f); fputc(0, f);
	// Write size data in big-endian 16-bit words
	fputc((image_extended_width & 0xFF00) >> 8, f);		// Extended width.
	fputc((image_extended_width & 0xFB), f);
	fputc((image_extended_height & 0xFF00) >> 8, f);	// Extended height.
	fputc((image_extended_height & 0xFB), f);
	fputc((image_width & 0xFF00) >> 8, f);		// Original width.
	fputc((image_width & 0xFF), f);
	fputc((image_height & 0xFF00) >> 8, f);		// Original height.
	fputc((image_height & 0xFF), f);
	int n = (image_extended_height / 4) * (image_extended_width / 4);
	fwrite(compressed_image, 1, n * 8, f);
	fclose(f);
}

static void write_dds_file() {
	printf("Writing .dds file.\n");
	int n = (image_extended_height / 4) * (image_extended_width / 4);
	FILE *f = fopen(output_filename, "wb");
	if (f == NULL) {
		printf("Error opening output file.\n");
		exit(1);
	}
	fputc('D', f); fputc('D', f); fputc('S', f); fputc(' ', f);
	unsigned char header[124];
	memset(header, 0, 124);
	*(unsigned int *)&header[0] = 124;
	*(unsigned int *)&header[4] = 0x1007;	// Flags
	*(unsigned int *)&header[8] = image_height;
	*(unsigned int *)&header[12] = image_width;
	*(unsigned int *)&header[16] = n * 8;	// Linear size.
	*(unsigned int *)&header[24] = 1;	// Mipmap count.
	*(unsigned int *)&header[72] = 32;
	*(unsigned int *)&header[76] = 0x4;	// Pixel format flags.
	header[80] = 'D';			// Pixel format.
	header[81] = 'X';
	header[82] = 'T';
	header[83] = '1';
	*(unsigned int *)&header[84] = 0x1000;	// Caps.
	fwrite(header, 1, 124, f);
	fwrite(compressed_image, 1, n * 8, f);
	fclose(f);
}

unsigned char ktx_id[12] = { 0xAB, 0x4B, 0x54, 0x58, 0x20, 0x31, 0x31, 0xBB, 0x0D, 0x0A, 0x1A, 0x0A };
char ktx_orientation_key[24] = { 'K', 'T', 'X', 'o', 'r', 'i', 'e', 'n', 't', 'a', 't', 'i', 'o', 'n', 0,
	'S', '=', 'r', ',', 'T', '=', 'd', 0, 0 };	// Includes one byte of padding.

static void write_ktx_file() {
	printf("Writing .ktx file.\n");
	int n = (image_extended_height / 4) * (image_extended_width / 4);
	FILE *f = fopen(output_filename, "wb");
	if (f == NULL) {
		printf("Error opening output file.\n");
		exit(1);
	}
	unsigned int header[16];
	memset(header, 0, 64);
	memcpy(header, ktx_id, 12);	// Set id.
	header[3] = 0x04030201;
	header[4] = 0;			// glType
	header[5] = 1;			// glTypeSize
	header[6] = 0;			// glFormat
	if (texture_format == 0)	
		header[7] = 0x8D64;		// GL_ETC1_RGB8_OES
	else
		header[7] = 0x9274;		// COMPRESSED_RGB8_ETC2
	header[9] = image_width;
	header[10] = image_height;
	header[11] = 0;
	header[13] = 1;		// Number of faces.
	header[14] = 1;		// Mipmap levels.
	int data[1];
	if (ktx_no_key_and_value) {
		header[15] = 0;
		fwrite(header, 1, 64, f);
	}
	else {
		header[15] = 28;	// Key value data bytes.
		fwrite(header, 1, 64, f);
		data[0] = 27;		// Key and value size.
		fwrite(data, 1, 4, f);
		fwrite(ktx_orientation_key, 1, 24, f);
	}
	data[0] = n * 8;	// Image size.
	fwrite(data, 1, 4, f);
	fwrite(compressed_image, 1, n * 8, f);
	fclose(f);
}

// Calculate the average difference per pixel of the whole image.

void calculate_difference() {
	double diff = 0;
	for (int y = 0; y < image_extended_height; y += 4)
		for (int x = 0; x < image_extended_width - 3; x+= 4) {
			sub_x_offset = x;
			sub_y_offset = y;
			double f = block4x4_calculate_fitness(NULL,
				(unsigned char *)&compressed_image[((y / 4) * (image_extended_width / 4) + x / 4) * 2]);
			if (f == 0) {
				printf("Error --unexpected fitness of 0 in calculate_difference.\n");
				exit(1);
			}
			diff += 1.0 / f;
		}
	int n = (image_extended_height / 4) * (image_extended_width / 4);
	double rmse = sqrt(diff / (n * 16));
	printf("Root-mean-square error per pixel: %lf  ", rmse);
	double psnr = 10.0 * log(255 * 255 / (diff / (3.0 * n * 16))) / log(10.0);
	printf("PSNR: %lf\n", psnr);
}

static void run_with_single_population() {
	printf("Running single GA for each pixel block.\n");
	FgenPopulation *pop = fgen_create(
		128,		// Population size.
		64,		// Number of bits.
		1,		// Data element size.
		block4x4_generation_callback,
		block4x4_calculate_fitness,
		fgen_seed_random,
		fgen_mutation_per_bit_fast,
		fgen_crossover_uniform_per_bit
		);
	fgen_set_parameters(
		pop,
		FGEN_ELITIST_SUS,
		FGEN_SUBTRACT_MIN_FITNESS,
		0.6,		// Crossover prob.
		0.01,		// Mutation prob. per bit
		0		// Macro-mutation prob.
		);
	fgen_set_generation_callback_interval(pop, nu_generations_per_subdivision);
	fgen_random_seed_with_timer(fgen_get_rng(pop));
//	fgen_random_seed_rng(fgen_get_rng(pop), 0);
	pop->user_data = malloc(sizeof(PopUserData));
	if (texture_format == 0)
		((PopUserData *)pop->user_data)->etc_modes_allowed = ETC_MODE_ALLOWED_ALL;
	else
	if (texture_format == 2)
		((PopUserData *)pop->user_data)->etc_modes_allowed = ETC2_MODE_ALLOWED_ALL;
	for (int y = 0; y < image_extended_height; y += 4)
		for (int x = 0; x < image_extended_width; x+= 4) {
			sub_x_offset = x;
			sub_y_offset = y;
			// Threading a single population doesn't work well in this case because the of the relatively quick
			// fitness function, which causes too much overhead caused by threading.
			if (use_threading)
				fgen_run_threaded(pop, nu_generations_per_subdivision);
			else
				fgen_run(pop, nu_generations_per_subdivision);
			report_solution(fgen_best_individual_of_population(pop));
		}
	free(pop->user_data);
	fgen_destroy(pop);
}

// The following function runs nu_islands concurrent GAs for the same pixel block. This doesn't improve speed but improves
// the accuracy (without increasing the number of generations) because of the wider range of solutions found.

static void run_with_archipelago() {
	FgenPopulation *pops[nu_islands];
	printf("Running GA archipelago of size %d for each pixel block.\n", nu_islands);
	for (int i = 0; i < nu_islands; i++) {
		pops[i] = fgen_create(
			128,		// Population size.
			64,		// Number of bits.
			1,		// Data element size.
			block4x4_generation_callback,
			block4x4_calculate_fitness,
			fgen_seed_random,
			fgen_mutation_per_bit_fast,
			fgen_crossover_uniform_per_bit
			);
		fgen_set_parameters(
			pops[i],
			FGEN_ELITIST_SUS,
			FGEN_SUBTRACT_MIN_FITNESS,
			0.6,		// Crossover prob.
			0.01,		// Mutation prob. per bit
			0		// Macro-mutation prob.
			);
		fgen_set_generation_callback_interval(pops[i], nu_generations_per_subdivision);
		fgen_set_migration_interval(pops[i], 0);	// No migration.
		fgen_set_migration_probability(pops[i], 0.01);
		pops[i]->user_data = malloc(sizeof(PopUserData));
		if (texture_format == 0) {
			if (modal && nu_islands >= 4)  {
				switch (i & 3) {
				case 0 :
				case 1 :
				case 2 :
					((PopUserData *)pops[i]->user_data)->etc_modes_allowed = ETC_MODE_ALLOWED_INDIVIDUAL;
						break;
				case 3 :
					((PopUserData *)pops[i]->user_data)->etc_modes_allowed = ETC_MODE_ALLOWED_DIFFERENTIAL;
						break;
				}
			}
			else
				((PopUserData *)pops[i]->user_data)->etc_modes_allowed = ETC_MODE_ALLOWED_ALL;
		}
		else
		if (texture_format == 2) {
			if (modal && nu_islands >= 8) {
				switch (i & 7) {
				case 0 :
				case 1 :
				case 2 :
					((PopUserData *)pops[i]->user_data)->etc_modes_allowed = ETC_MODE_ALLOWED_INDIVIDUAL;
					break;
				case 3 :
				case 4 :
					((PopUserData *)pops[i]->user_data)->etc_modes_allowed = ETC_MODE_ALLOWED_DIFFERENTIAL;
					break;
				case 5 :
					((PopUserData *)pops[i]->user_data)->etc_modes_allowed = ETC2_MODE_ALLOWED_T;
					break;
				case 6 :
					((PopUserData *)pops[i]->user_data)->etc_modes_allowed = ETC2_MODE_ALLOWED_H;
					break;
				case 7 :
					((PopUserData *)pops[i]->user_data)->etc_modes_allowed = ETC2_MODE_ALLOWED_PLANAR;
					break;
				}
			}
			else
				((PopUserData *)pops[i]->user_data)->etc_modes_allowed = ETC2_MODE_ALLOWED_ALL;
		}
	}
	fgen_random_seed_with_timer(fgen_get_rng(pops[0]));
//	fgen_random_seed_rng(fgen_get_rng(pops[0]), 0);

	for (int y = 0; y < image_extended_height; y += 4)
		for (int x = 0; x < image_extended_width; x+= 4) {
			sub_x_offset = x;
			sub_y_offset = y;
			if (use_threading)
				fgen_run_archipelago_threaded(nu_islands, pops, nu_generations_per_subdivision);
			else
				fgen_run_archipelago(nu_islands, pops, nu_generations_per_subdivision);
			if (verbose && (texture_format == 0 || texture_format == 2)) {
				for (int i = 0; i < nu_islands; i++) {
					printf("Block %d: ", (y / 4) * (image_extended_width / 4) + (x / 4));
					printf("Modes: ");
					int modes_allowed =  ((PopUserData *)pops[i]->user_data)->etc_modes_allowed;
					if (modes_allowed & ETC_MODE_ALLOWED_INDIVIDUAL)
						printf("I");
					if (modes_allowed & ETC_MODE_ALLOWED_DIFFERENTIAL)
						printf("D");
					if (modes_allowed & ETC2_MODE_ALLOWED_T)
						printf("T");
					if (modes_allowed & ETC2_MODE_ALLOWED_H)
						printf("H");
					if (modes_allowed & ETC2_MODE_ALLOWED_PLANAR)
						printf("P");
					FgenIndividual *best = fgen_best_individual_of_population(pops[i]);
					printf(" RMSE per pixel: %lf\n", sqrt((1.0 /
						block4x4_calculate_fitness(NULL, best->bitstring)) / 16));
				}
			}
			FgenIndividual *best = fgen_best_individual_of_archipelago(nu_islands, pops);
			if (verbose)
				printf("Block %d: Combined: RMSE per pixel: %lf\n", (y / 4 ) * (image_extended_width / 4) + (x / 4),
					sqrt((1.0 / block4x4_calculate_fitness(NULL, best->bitstring)) / 16));
			report_solution(best);
		}
	for (int i = 0; i < nu_islands; i++) {
		free(pops[i]->user_data);
		fgen_destroy(pops[i]);
	}
}


void start_ga_block4x4() {
	clear_result_image();

	sub_image_buffer = malloc(128);
	compressed_image = malloc((image_extended_height / 4) * (image_extended_width / 4) * 8);

	switch (texture_format) {
	case 0 :
		printf("Using ETC1 compression.\n");
		break;
	case 1 :
		printf("Using DXT1 compression.\n");
		break;
	case 2 :
		printf("Using ETC2 (RGB8) compression.\n");
		break;
	}

	if (nu_islands > 1)
		run_with_archipelago();
	else
		run_with_single_population();

	update_area2();
	calculate_difference();

	if (output_filename != NULL) {
		// Write the compressed image.
		switch (texture_output_format) {
		case 0 :
			write_pkm_file();
			break;
		case 1 :
			write_dds_file();
			break;
		case 2 :
			write_ktx_file();
			break;
		}
		printf("Succesfully wrote compressed texture to the file %s.\n", output_filename);
	}
}

static void read_pkm_file() {
	printf("Reading .pkm file %s.\n", output_filename);
	FILE *f = fopen(output_filename, "r");
	if (f == NULL) {
		printf("Error opening file.\n");
		exit(1);
	}
	unsigned char header[16];
	fread(header, 1, 16, f);
	if (header[0] != 'P' || header[1] != 'K' || header[2] != 'M' || header[3] != ' ') {
		printf("Error -- couldn't find PKM signature.\n");
		exit(1);
	}
	int texture_type = ((int)header[6] << 8) | header[7];
	if (texture_type != 0 && texture_type != 1) {
		printf("Error -- unsupported format (only ETC1 and ETC2_RGB supported).\n");
		exit(1);
	}
	if (texture_type == 1 && texture_format == 0) {
		printf("Switching selected texture format from ETC1 to ETC2.\n");
		texture_format = 2;
	}
	int ext_width = ((int)header[8] << 8) | header[9];
	int ext_height = ((int)header[10] << 8) | header[11];
	if (ext_width != image_extended_width || ext_height != image_extended_height) {
		printf("Error -- size mismatch between files.\n");
		exit(1);
	}
	int n = (image_height / 4) * (image_width / 4);
	if (fread(compressed_image, 1, n * 8, f) < n * 8) {
		printf("Error reading file.\n");
		exit (1);
	}
	fclose(f);
}


static void read_dds_file() {
	printf("Reading .dds file %s.\n", output_filename);
	FILE *f = fopen(output_filename, "r");
	if (f == NULL) {
		printf("Error opening file.\n");
		exit(1);
	}
	char id[4];
	fread(id, 1, 4, f);
	if (id[0] != 'D' || id[1] != 'D' || id[2] != 'S' || id[3] != ' ') {
		printf("Error -- couldn't find DDS signature.\n");
		exit(1);
	}
	unsigned char header[124];
	fread(header, 1, 124, f);
	if (*(unsigned int *)&header[8] != image_height ||
	*(unsigned int *)&header[12] != image_width) {
		printf("Error -- size mismatch between files.\n");
		exit(1);
	}
	if (header[80] != 'D' || header[81] != 'X' || header[82] != 'T' || header[83] != '1') {
		printf("Error -- unsupported format (only DXT1 supported).\n");
		exit(1);
	}
	// Only read the first mip-map level.
	int n = (image_extended_height / 4) * (image_extended_width / 4);
	int r = fread(compressed_image, 1, n * 8, f);
	if (r < n * 8) {
		printf("Error reading file.\n");
		exit(1);
	}
	fclose(f);
}

static void read_ktx_file() {
	printf("Reading .ktx file %s.\n", output_filename);
	FILE *f = fopen(output_filename, "r");
	if (f == NULL) {
		printf("Error opening file.\n");
		exit(1);
	}
	int header[16];
	fread(header, 1, 64, f);
	if (memcmp(header, ktx_id, 12) != 0) {
		printf("Error -- couldn't find KTX signature.\n");
		exit(1);
	}
	int big_endian = 0;
	if (header[3] == 0x01020304) {
		// Big endian .ktx file.
		big_endian = 1;
		for (int i = 3; i < 16; i++) {
			unsigned char *b = (unsigned char *)&header[i];
			unsigned char temp = b[0];
			b[0] = b[3];
			b[3] = temp;
			temp = b[1];
			b[1] = b[2];
			b[2] = temp;
		}
	}
	if (header[7] != 0x8D64 && header[7] != 0x9274) {	// GL_ETC1_RGB8_OES
		printf("Error -- unsupported format (only ETC1 and RGB8 ETC2 supported).\n");
		exit(1);
	}
	if (header[7] == 0x9274 && texture_format == 0) {
		printf("Switching selected texture format from ETC1 to ETC2.\n");
		texture_format = 2;
	}
	if (header[9] != image_width || header[10] != image_height) {
		printf("Error -- size mismatch between files.\n");
		exit(1);
	}
	if (header[15] > 0) {
		// Skip metadata.
		unsigned char *metadata = malloc(header[15]);
		if (fread(metadata, 1, header[15], f) < header[15]) {
			printf("Error reading metadata.\n");
			exit(1);
		}
		free(metadata);
	}
	int n = (image_extended_height / 4) * (image_extended_width / 4);
	unsigned char image_size[4];
	fread(image_size, 1, 4, f);
	if (big_endian) {
		unsigned char temp = image_size[0];
		image_size[0] = image_size[3];
		image_size[3] = temp;
		temp = image_size[1];
		image_size[1] = image_size[2];
		image_size[2] = temp;
	}
	if (*(int *)&image_size[0] != n * 8) {
		printf("Error -- image size field does not match (%d vs %d).\n", *(int *)&image_size[0], n * 8);
		exit(1);
	}
	if (fread(compressed_image, 1, n * 8, f) < n * 8) {
		printf("Error reading file.\n");
		exit (1);
	}
	fclose(f);
}

void compare_files() {
	compressed_image = malloc((image_extended_height / 4) * (image_extended_width / 4) * 8);
	switch (texture_output_format) {
	case 0 :
		read_pkm_file();
		break;
	case 1 :
		read_dds_file();
		break;
	case 2 :
		read_ktx_file();
		break;
	}
	calculate_difference();
}

