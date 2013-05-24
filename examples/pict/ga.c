/*
    ga.c -- genetic algorithm interface of pict, a graphical example using fgen.

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
#include <stdio.h>
#include <string.h>
#include <gtk/gtk.h>
#include <fgen.h>

#include "pict.h"

const int NU_PARAMETERS_PER_COMPONENT = 11;
int nu_components = 2;	// Configurable.
int nu_generations_per_subdivision = 0;	// Configurable.
static int nu_params;
static Ffit *fit;
static unsigned int *sub_image_buffer;
static int sub_image_width;
static int sub_image_height;
static int sub_x_offset, sub_y_offset;
static int max_generation;

/*
	Grid mode:
	Background:
	0	red
	1	green
	2	blue
	3	type [0, 0,5[ = constant-filled
		     [0.5, 0.75[ = horizontal gradient-filled
		     [0.75, 1.0[ = vertical gradient-filled
	4	second red component for gradient
	5	second green component for gradient
	6	second blue component for gradient
	<components>
	Parameters per component:
	0	x-coordinate in pixels
	1	y-coordinate in pixels
	2	width in pixels
	3	height in pixels
	4	red component [0, 1.0[
	5	green component
	6	blue component
	7	type [0, 0,5[ = constant-filled
		     [0.5, 0.75[ = horizontal gradient-filled
		     [0.75, 1.0[ = vertical gradient-filled
	8	second red component for gradient
	9	second green component for gradient
	10	second blue component for gradient

	Grid mode:
	<grid colors>
	Parameters per color:
	0	red component
	1	green component
	2	blue component
*/

static void draw_rectangle(unsigned int *image_buffer, int x, int y, int w, int h, int r, int g, int b) {
	if (w == 0 || h == 0)
		return;
	unsigned int *pixels = image_buffer;
	unsigned int value = r + (g << 8) + (b << 16);
	int ye = y + h;
	if (ye > sub_image_height)
		ye = sub_image_height;
	int xe = x + w;
	if (xe > sub_image_width)
		xe = sub_image_width;
	unsigned int *pix = &pixels[y * sub_image_width];
	for (int yy = y; yy < ye; yy++) {
		int xx = x;
		while (xe - xx >= 8) {
			pix[xx] = value;
			pix[xx + 1] = value;
			pix[xx + 2] = value;
			pix[xx + 3] = value;
			pix[xx + 4] = value;
			pix[xx + 5] = value;
			pix[xx + 6] = value;
			pix[xx + 7] = value;
			xx += 8;
		}
		for (; xx < xe; xx++)
			pix[xx] = value;
		pix = &pix[sub_image_width];
	}
}

static void draw_rectangle_vertical_gradient(unsigned int *image_buffer, int x, int y, int w, int h,
int r, int g, int b, int r2, int g2, int b2) {
	if (h == 0  || w == 0)
		return;
	unsigned int *pixels = image_buffer;
	float r_step = (float)(r2 - r) / h;
	float g_step = (float)(g2 - g) / h;
	float b_step = (float)(b2 - b) / h;
	int ye = y + h;
	if (ye > sub_image_height)
		ye = sub_image_height;
	int xe = x + w;
	if (xe > sub_image_width)
		xe = sub_image_width;
	unsigned int *pix = &pixels[y * sub_image_width];
	for (int yy = y; yy < ye; yy++) {
		int red = r + r_step * (yy - y);
		int green = g + g_step * (yy - y);
		int blue = b + b_step * (yy - y);
		unsigned int value = red + (green << 8) + (blue << 16);
		int xx = x;
		while (xe - xx >= 8) {
			pix[xx] = value;
			pix[xx + 1] = value;
			pix[xx + 2] = value;
			pix[xx + 3] = value;
			pix[xx + 4] = value;
			pix[xx + 5] = value;
			pix[xx + 6] = value;
			pix[xx + 7] = value;
			xx += 8;
		}
		for (; xx < xe; xx++)
			pix[xx] = value;
		pix = &pix[sub_image_width];
	}
}

static void draw_rectangle_horizontal_gradient(unsigned int *image_buffer, int x, int y, int w, int h,
int r, int g, int b, int r2, int g2, int b2) {
	if (w == 0 || h == 0)
		return;
	unsigned int *pixels = image_buffer;
	float r_step = (float)(r2 - r) / w;
	float g_step = (float)(g2 - g) / w;
	float b_step = (float)(b2 - b) / w;
	int xe = x + w;
	if (xe > sub_image_width)
		xe = sub_image_width;
	int ye = y + h;
	if (ye > sub_image_height)
		ye = sub_image_height;
	for (int xx = x; xx < xe; xx++) {
		int red = r + r_step * (xx - x);
		int green = g + g_step * (xx - x);
		int blue = b + b_step * (xx - x);
		unsigned int value = red + (green << 8) + (blue << 16);
		unsigned int *pix = &pixels[y * sub_image_width];
		for (int yy = y; yy < ye; yy++) {
			pix[xx] = value;
			pix = &pix[sub_image_width];
		}
	}
}


int compare_area(const void *e1, const void *e2) {
	double *component1 = *(double **)e1;
	double *component2 = *(double **)e2;
	// Calculate the clipped rectangle areas.
	int x1 = component1[0];
	int y1 = component1[1];
	int w1 = component1[2];
	int h1 = component1[3];
	if (x1 + w1 > sub_image_width)
		w1 = sub_image_width - x1;
	if (y1 + h1 > sub_image_height)
		h1 = sub_image_height - y1;
	int x2 = component2[0];
	int y2 = component2[1];
	int w2 = component2[2];
	int h2 = component2[3];
	if (x2 + w2 > sub_image_width)
		w2 = sub_image_width - x2;
	if (y2 + h2 > sub_image_height)
		h2 = sub_image_height - y2;
	int area1 = w1 * h1;
	int area2 = w2 * h2;
	if (area1 > area2)
		return - 1;
	if (area1 < area2);
		return 1;
	return 0;
}

static void draw_picture(unsigned int *image_buffer, const double *param) {
	// Clear picture with background color.
	if (param[3] >= 0.5) {
		// Gradient.
		int r2 = param[4] * 255.5;
		int g2 = param[5] * 255.5;
		int b2 = param[6] * 255.5;
		if (param[3] >= 0.75)
			draw_rectangle_vertical_gradient(image_buffer, 0, 0, sub_image_width, sub_image_height,
				param[0] * 255.5, param[1] * 255.5, param[2] * 255.5, r2, g2, b2);
		else
			draw_rectangle_horizontal_gradient(image_buffer, 0, 0, sub_image_width, sub_image_height,
				param[0] * 255.5, param[1] * 255.5, param[2] * 255.5, r2, g2, b2);
	}
	else
		draw_rectangle(image_buffer, 0, 0, sub_image_width, sub_image_height, param[0] * 255.5,
			param[1] * 255.5, param[2] * 255.5);
	// Sort the rectangles on decreasing area.
	double **sorted = alloca(sizeof(double *) * nu_components);
	for (int i = 0; i < nu_components; i++)
		sorted[i] = &param[i * NU_PARAMETERS_PER_COMPONENT + 7];
	qsort(sorted, nu_components, sizeof(double *), compare_area);
	// Draw all component filled rectangles.
	for (int i = 0; i < nu_components; i++) {
		int x = sorted[i][0];
		int y = sorted[i][1];
		int w = sorted[i][2];
		int h = sorted[i][3];
		int r = sorted[i][4] * 255.5;
		int g = sorted[i][5] * 255.5;
		int b = sorted[i][6] * 255.5;
		if (sorted[i][7] >= 0.5) {
			// Gradient.
			int r2 = sorted[i][8] * 255.5;
			int g2 = sorted[i][9] * 255.5;
			int b2 = sorted[i][10] * 255.5;
			if (sorted[i][7] >= 0.75)
				draw_rectangle_vertical_gradient(image_buffer, x, y, w, h, r, g, b, r2, g2, b2);
			else
				draw_rectangle_horizontal_gradient(image_buffer, x, y, w, h, r, g, b, r2, g2, b2);
		}
		else {
			// Constant-filled rectangle.
			draw_rectangle(image_buffer, x, y, w, h, r, g, b);
		}
	}
}


// Grid mode.

static void draw_picture_grid(unsigned int *image_buffer, const double *param) {
	int i = 0;
	for (int y = 0; y < grid_size; y++)
		for (int x = 0; x < grid_size; x++) {
			int xx = sub_image_width * x / grid_size;
			int xx_next = sub_image_width * (x + 1) / grid_size;
			int yy = sub_image_height * y / grid_size;
			int yy_next = sub_image_height * (y + 1) / grid_size;
			int w = xx_next - xx;
			int h = yy_next - yy;
			int r = param[i] * 255.5;
			int g = param[i + 1] * 255.5;
			int b = param[i + 2] * 255.5;	
			draw_rectangle(image_buffer, xx, yy, w, h, r, g, b);
			i += 3;
		}
}

static void print_params(const double *param) {
	// Sort the rectangle on decreasing area.
	double **sorted = malloc(sizeof(double *) * nu_components);
	for (int i = 0; i < nu_components; i++)
		sorted[i] = &param[i * NU_PARAMETERS_PER_COMPONENT + 7];
	qsort(sorted, nu_components, sizeof(double *), compare_area);
	for (int i = 0; i < nu_components; i++) {
		// Calculate the clipped rectangle areas.
		int x1 = sorted[i][0];
		int y1 = sorted[i][1];
		int w1 = sorted[i][2];
		int h1 = sorted[i][3];
		if (x1 + w1 > sub_image_width)
			w1 = sub_image_width - x1;
		if (y1 + h1 > sub_image_height)
			h1 = sub_image_height - y1;
		int area = w1 * h1;
		printf("Rect (%d, %d, %d, %d) area = %d.\n", x1, y1, w1, h1, area);
	}
	free(sorted);
}

static double calculate_error(const Ffit *fit, const double *param) {
	if (grid_size > 0)
		draw_picture_grid(sub_image_buffer, param);
	else
		draw_picture(sub_image_buffer, param);
	unsigned int *pix = sub_image_buffer;
	guchar *pixbuf_pix = image_pixels + sub_y_offset * image_rowstride + sub_x_offset * 3;
	int error = 0;
	for (int y = 0; y < sub_image_height; y++) {
		for (int x = 0; x < sub_image_width; x++) {
			unsigned char *pixb = (unsigned char *)&pix[x];
			unsigned int r = pixb[0];
			unsigned int g = pixb[1];
			unsigned int b = pixb[2];
			unsigned int r_source = pixbuf_pix[x * 3];
			unsigned int g_source = pixbuf_pix[x * 3 + 1];
			unsigned int b_source = pixbuf_pix[x * 3 + 2];
			error += abs(r - r_source);
			error += abs(g - g_source);
			error += abs(b - b_source);
		}
		pix = &pix[sub_image_width];
		pixbuf_pix = &pixbuf_pix[image_rowstride];
	}
	return error;
}

static double calculate_error_thread(const Ffit *fit, const double *param) {
	unsigned int *image_buffer = alloca(sub_image_height * sub_image_width * 4);
	if (grid_size > 0)
		draw_picture_grid(image_buffer, param);
	else
		draw_picture(image_buffer, param);
	unsigned int *pix = image_buffer;
	guchar *pixbuf_pix = image_pixels + sub_y_offset * image_rowstride + sub_x_offset * 3;
	int error = 0;
	for (int y = 0; y < sub_image_height; y++) {
		for (int x = 0; x < sub_image_width; x++) {
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
		pix = &pix[sub_image_width];
		pixbuf_pix = &pixbuf_pix[image_rowstride];
	}
	return error;
}

void copy_buffer_to_result_pixbuf(unsigned int *buffer, int sub_image_width, int sub_image_height, int sub_x_offset,
int sub_y_offset) {
	guchar *pixels_result = gdk_pixbuf_get_pixels(result_pixbuf);
	for (int y = 0; y < sub_image_height; y++)
		for (int x = 0; x < sub_image_width; x++) {
			unsigned int value = buffer[y * sub_image_width + x];
			pixels_result[(y + sub_y_offset) * image_rowstride + (x + sub_x_offset) * 3] =
				value & 0xFF;
			pixels_result[(y + sub_y_offset) * image_rowstride + (x + sub_x_offset) * 3 + 1] =
				(value & 0xFF00) >> 8;
			pixels_result[(y + sub_y_offset) * image_rowstride + (x + sub_x_offset) * 3 + 2] =
				(value & 0xFF0000) >> 16;
		}
}

static void generation_callback(Ffit *fit, int generation, const double *best_param, double best_error) {
	while (gtk_events_pending())
		gtk_main_iteration();
	if (generation == 0)
		return;

//	printf("Generation %d, best error = %lf.\n", generation, best_error);
	int n = nu_subdivisions;
	if (grid_size > 0)
		n = 1;
	if (n == 1 && (generation % 2) != 0)
		return;
	if (n == 2 && (generation % 5) != 0)
		return;
	if (n == 4 && (generation % 25) != 0)
		return;
	if (n > 4 && (generation % 100 != 0 && generation < max_generation))
		return;
	if (grid_size > 0)
		draw_picture_grid(sub_image_buffer, best_param);
	else
		draw_picture(sub_image_buffer, best_param);
//	print_params(best_param);
	copy_buffer_to_result_pixbuf(sub_image_buffer, sub_image_width, sub_image_height, sub_x_offset,
		sub_y_offset);
	update_area2();
	while (gtk_events_pending())
		gtk_main_iteration();
	if (generation >= max_generation)
		ffit_signal_stop(fit);
}

void start_ga(int n) {
	nu_params = 7 + nu_components * NU_PARAMETERS_PER_COMPONENT;
	clear_result_image();

	for (int y = 0; y < n; y++)
		for (int x = 0; x < n; x++) {
			sub_x_offset = image_width * x / n;
			sub_y_offset = image_height * y / n;
			int sub_x_offset_next = image_width * (x + 1) / n;
			int sub_y_offset_next = image_height * (y + 1) / n;
			sub_image_width = sub_x_offset_next - sub_x_offset;
			sub_image_height = sub_y_offset_next - sub_y_offset;
			sub_image_buffer = malloc(4 * sub_image_width * sub_image_height);

			if (use_threading) {
				fit = ffit_create(nu_params, generation_callback, calculate_error_thread);
				ffit_set_threading(fit, FFIT_THREADING_ENABLED);
			}
			else
				fit = ffit_create(nu_params, generation_callback, calculate_error);
			int j = 0;
			// Background color.
			ffit_set_parameter_range_and_mapping(fit, j, 0, 1.0, FFIT_MAPPING_LINEAR);
			ffit_set_parameter_range_and_mapping(fit, j + 1, 0, 1.0, FFIT_MAPPING_LINEAR);
			ffit_set_parameter_range_and_mapping(fit, j + 2, 0, 1.0, FFIT_MAPPING_LINEAR);
			ffit_set_parameter_range_and_mapping(fit, j + 3, 0, 1.0, FFIT_MAPPING_LINEAR);
			ffit_set_parameter_range_and_mapping(fit, j + 4, 0, 1.0, FFIT_MAPPING_LINEAR);
			ffit_set_parameter_range_and_mapping(fit, j + 5, 0, 1.0, FFIT_MAPPING_LINEAR);
			ffit_set_parameter_range_and_mapping(fit, j + 6, 0, 1.0, FFIT_MAPPING_LINEAR);
			j += 7;
			// Components.
			for (int i = 0; i < nu_components; i++) {
				ffit_set_parameter_range_and_mapping(fit, j, 0, sub_image_width,
					FFIT_MAPPING_LINEAR);
				ffit_set_parameter_range_and_mapping(fit, j + 1, 0, sub_image_height,
					FFIT_MAPPING_LINEAR);
				ffit_set_parameter_range_and_mapping(fit, j + 2, 0, sub_image_width + 0.5,
					FFIT_MAPPING_SQUARE);
				ffit_set_parameter_range_and_mapping(fit, j + 3, 0, sub_image_height + 0.5,
					FFIT_MAPPING_SQUARE);
				ffit_set_parameter_range_and_mapping(fit, j + 4, 0, 1.0, FFIT_MAPPING_LINEAR);
				ffit_set_parameter_range_and_mapping(fit, j + 5, 0, 1.0, FFIT_MAPPING_LINEAR);
				ffit_set_parameter_range_and_mapping(fit, j + 6, 0, 1.0, FFIT_MAPPING_LINEAR);
				ffit_set_parameter_range_and_mapping(fit, j + 7, 0, 1.0, FFIT_MAPPING_LINEAR);
				ffit_set_parameter_range_and_mapping(fit, j + 8, 0, 1.0, FFIT_MAPPING_LINEAR);
				ffit_set_parameter_range_and_mapping(fit, j + 9, 0, 1.0, FFIT_MAPPING_LINEAR);
				ffit_set_parameter_range_and_mapping(fit, j + 10, 0, 1.0, FFIT_MAPPING_LINEAR);
				j += NU_PARAMETERS_PER_COMPONENT;
			}
			max_generation = nu_generations_per_subdivision;
			ffit_run_fgen(fit, 128, 16, FGEN_ELITIST_SUS, fgen_crossover_uniform_per_element, 0.8,
				0.015, 0.05);
			ffit_destroy(fit);
			free(sub_image_buffer);
		}
}

void start_ga_grid(int n) {
	nu_params = grid_size * grid_size * 3;
	clear_result_image();

	for (int y = 0; y < n; y++)
		for (int x = 0; x < n; x++) {
			sub_x_offset = image_width * x / n;
			sub_y_offset = image_height * y / n;
			int sub_x_offset_next = image_width * (x + 1) / n;
			int sub_y_offset_next = image_height * (y + 1) / n;
			sub_image_width = sub_x_offset_next - sub_x_offset;
			sub_image_height = sub_y_offset_next - sub_y_offset;
//			sub_pixels = gdk_pixbuf_get_pixels(scratch_pixbuf) + sub_y_offset * image_rowstride +
//				 sub_x_offset * 3;
			sub_image_buffer = malloc(4 * sub_image_width * sub_image_height);

			if (use_threading) {
				fit = ffit_create(nu_params, generation_callback, calculate_error_thread);
				ffit_set_threading(fit, FFIT_THREADING_ENABLED);
			}
			else
				fit = ffit_create(nu_params, generation_callback, calculate_error);
			int j = 0;
			int k = grid_size * grid_size;
			for (int i = 0; i < k; i++) {
				ffit_set_parameter_range_and_mapping(fit, j, 0, 1.0, FFIT_MAPPING_LINEAR);
				ffit_set_parameter_range_and_mapping(fit, j + 1, 0, 1.0, FFIT_MAPPING_LINEAR);
				ffit_set_parameter_range_and_mapping(fit, j + 2, 0, 1.0, FFIT_MAPPING_LINEAR);
				j += 3;
			}
			max_generation = nu_generations_per_subdivision;
			ffit_run_fgen(fit, 128, 16, FGEN_ELITIST_SUS, fgen_crossover_uniform_per_element, 0.8,
				0.015, 0.05);
			ffit_destroy(fit);
			free(sub_image_buffer);
		}
}


