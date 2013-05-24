/*
    main.c -- main file of pict, a graphical example using fgen.

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

// Define _GNU_SOURCE for strcasestr.
#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <gtk/gtk.h>
#include "pict.h"

#ifdef _WIN32

char *strcasestr(char *haystack, char *needle) {
        char *p, *startn = NULL, *np = NULL;

        for (p = haystack; *p != '\0'; p++) {
                if (np) {
                        if (toupper(*p) == toupper(*np)) {
                                np++;
                                if (*np == '\0')
                                        return startn;
                        } else
                                np = NULL;
                } else if (toupper(*p) == toupper(*needle)) {
                        np = needle + 1;
                        startn = p;
                }
        }
        return NULL;
}

#endif


GdkPixbuf *image_pixbuf, *result_pixbuf, *magnified_pixbuf;
int image_width, image_height, image_rowstride;
int image_extended_width, image_extended_height;
guchar *image_pixels;
int nu_subdivisions = 1;
int grid_size = 0;
static int auto_quit = FALSE;
int use_threading = FALSE;
int block4x4_mode = FALSE;
int magnification = 1;
char *output_filename = NULL;
int texture_output_format = 0;
int texture_format = 0;
int ktx_no_key_and_value = 0;
int nu_islands = 0;
static int compare_files_flag = 0;
int verbose = 0;
int modal = 0;

void load_image_file(char *filename) {
	GdkPixbufFormat *format;
	int width, height;
	format = gdk_pixbuf_get_file_info(filename, &width, &height);
	if (format == NULL) {
		printf("Error - file %s is not a recognized image file.\n", filename);
		exit(1);
	}
	GError *error = NULL;
	GdkPixbuf *pixbuf = gdk_pixbuf_new_from_file(filename, &error);
	if (pixbuf == NULL) {
		printf("gdk_pixbuf_new_from_file returned error: %s.\n", error->message);
		exit(1);
	}
	int bits_per_sample = gdk_pixbuf_get_bits_per_sample(pixbuf);
	if (bits_per_sample != 8) {
		printf("Error -- expected 8 bits per sample.\n");
		exit(1);
	}
	image_width = gdk_pixbuf_get_width(pixbuf);
	image_height = gdk_pixbuf_get_height(pixbuf);
	image_pixels = gdk_pixbuf_get_pixels(pixbuf);
	image_rowstride = gdk_pixbuf_get_rowstride(pixbuf);
	image_extended_width = (image_width + 3) & ~3;
	image_extended_height = (image_height + 3) & ~3;
	if (gdk_pixbuf_get_has_alpha(pixbuf) == TRUE) {
		printf("Discarding alpha channel.\n");
		GdkPixbuf *new = gdk_pixbuf_new(GDK_COLORSPACE_RGB, FALSE, 8, image_width, image_height);
		guchar *new_pixels = gdk_pixbuf_get_pixels(new);
		int new_rowstride = gdk_pixbuf_get_rowstride(new);
		for (int y = 0; y < image_height; y++)
			for (int x = 0; x < image_width; x++) {
				memcpy(new_pixels + y * new_rowstride + x * 3,
					image_pixels + y * image_rowstride + x * 4, 3);
			}
		g_object_unref(G_OBJECT(pixbuf));
		pixbuf = new;
		image_pixels = new_pixels;
		image_rowstride = new_rowstride;
	}
	printf("Loaded image %s, width %d, height %d, bits per sample %d, rowstride %d.\n",
		filename, image_width, image_height, bits_per_sample, image_rowstride);
	image_pixbuf = pixbuf;
}

void create_result_pixbuf() {
	result_pixbuf = gdk_pixbuf_copy(image_pixbuf);
	if (magnification > 1) {
		magnified_pixbuf = gdk_pixbuf_new(GDK_COLORSPACE_RGB, FALSE, 8, gdk_pixbuf_get_width(image_pixbuf) *
			magnification, gdk_pixbuf_get_height(image_pixbuf) * magnification);
	}
}

gboolean start_ga_cb(gpointer data) {
	if (block4x4_mode)
		start_ga_block4x4();
	else
	if (grid_size > 0)
		start_ga_grid(nu_subdivisions);
	else
		start_ga(nu_subdivisions);
	if (auto_quit)
		gtk_main_quit();
	return FALSE;
}

int main(int argc, char **argv) {
	gtk_init( &argc, &argv );
	if (argc == 1) {
		printf("pict -- genetic algorithm picture decompositor using rectangles (constant\n"
			"or gradient-filled), plus texture compression.\n"
			"Usage: pict <options> <image filename>\n"
			"Options:\n"
			"-n <number>   Number of horizontal and vertical subdivisions to use. Default 1.\n"
			"-c <number>   Number of geometrical components per subdivision.\n"
			"              Ignored if -f is used. Default 2.\n"
			"-g <number>   Max generations to run per subdivision if n > 1. Default 100.\n"
			"-f <number>   Use a fixed grid of n x n rectangles per subdivision instead\n"
			"              of arbitrary ones.\n"
			"-q            Quit when finished.\n"
#if __WORDSIZE == 64
			"-t            Use threading.\n"
#endif
			"-a <number>   In texture compression, use an archipelago of <number> concurrent GAs\n"
			"              for each block. Note that this doesn't speed up the process but does\n"
			"              improve the result. Enabling -t recommended, optimal number depends on\n"
			"              CPU cores and cache. 8 or 16 is a good choice for modern CPUs.\n"
			"-e            4x4 pixel ETC1 texture compression mode. Supports option -g.\n"
			"-e2           4x4 pixel ETC2 RGB texture compression mode. Supports option -g.\n"
			"-d            4x4 pixel DXT1 texture compression mode. Supports option -g.\n"
			"-m <number>   Magnify picture by factor (integer).\n"
			"-o <filename> Output filename. When specified with option -e or -d, the\n"
			"              compressed image is written to this file. The filename must\n"
			"              have the .pkm or .ktx extension for ETC1 compression or .dds\n"
			"              for DXT1 compression or .ktx for ETC2 compression.\n"
			"-k            When writing a .ktx file, don't use key and value pairs because\n"
			"              certain applications (such as etcpack) can't handle them.\n"
			"-x            Only compare the image file with the compressed 'output' file\n"
			"              specified by option -o and report the difference. No file is written.\n"
			"-v            Information printed for each block.\n"
			"-l            When doing ETC compression with an archipelago with at least 4 (ETC1)\n"
			"              or 8 (ETC2) islands, run each island limited to a single compression mode.\n"
			"\nPress 'q' to quit when running.\n");
		exit(0);
	}

	int i = 1;
	for (; i + 1 < argc;) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			case 'n' :
				nu_subdivisions = atoi(argv[i + 1]);
				i += 2;
				break;
			case 'c' :
				nu_components = atoi(argv[i + 1]);
				i += 2;
				break;
			case 'g' :
				nu_generations_per_subdivision = atoi(argv[i + 1]);
				i += 2;
				break;
			case 'f' :
				grid_size = atoi(argv[i + 1]);
				i += 2;
				break;
			case 'q' :
				auto_quit = TRUE;
				i++;
				break;
#if __WORDSIZE == 64
			case 't' :
				use_threading = TRUE;	
				i++;
				break;
#endif
			case 'a' :
				nu_islands = atoi(argv[i + 1]);
				i += 2;
				break;
			case 'e' :
				block4x4_mode = TRUE;
				if (argv[i][2] == '2')
					texture_format = 2;	// ETC2
				else
					texture_format = 0;	// ETC1
				i++;
				break;
			case 'd' :
				block4x4_mode = TRUE;
				texture_format = 1;	// DXT1
				i++;
				break;
			case 'm' :
				magnification = atoi(argv[i + 1]);
				i += 2;
				break;
			case 'o' :
				output_filename = argv[i + 1];
				if (strcasestr(output_filename, ".pkm") != NULL)
					texture_output_format = 0;
				else
				if (strcasestr(output_filename, ".dds") != NULL)
					texture_output_format = 1;
				else
				if (strcasestr(output_filename, ".ktx") != NULL)
					texture_output_format = 2;
				else {
					printf("Error -- expected filename with extension .pkx or .ktx (for etc1 compression) "
						"or .dds or (for DXT1 compression).\n");
					exit(1);
				}
				i += 2;
				break;
			case 'k' :
				ktx_no_key_and_value = 1;
				i++;
				break;
			case 'x' :
				compare_files_flag  = 1;
				i++;
				break;
			case 'v' :
				verbose = 1;
				i++;
				break;
			case 'l' :
				modal = 1;
				i++;
				break;
			default :
				printf("Invalid option.\n");
				exit(1);
			}
			continue;
		}
		break;
	}

	if (output_filename != NULL && ((texture_format == 0 && texture_output_format != 0 && texture_output_format != 2) ||
	(texture_format == 1 && texture_output_format != 1) || (texture_format == 2 && texture_output_format != 2))) {
		printf("Invalid combination of texture type and output filename type.\n");
		exit(1);
	}

	if (nu_generations_per_subdivision == 0) {
		if (nu_subdivisions == 1 && !block4x4_mode)
			nu_generations_per_subdivision = INT_MAX;
		else
			nu_generations_per_subdivision = 100;
	}

	if (i >= argc) {
		printf("Must specify image filename.\n");
		exit(1);
	}

	FILE *f;
	f = fopen(argv[i], "rb");
	if (f == NULL) {
		printf("Error reading file %s.\n", argv[i]);
		exit(1);
	}
	fclose(f);

	load_image_file(argv[i]);
	if (output_filename != NULL && ((image_width & 3) != 0 || (image_height & 3) != 0)) {
		printf("Error -- expected width and height to be a multiple of four for texture compression.\n");
		exit(1);
	}
	// Avoid window that is too large.
	if (image_width * magnification > 1920 * 2 || image_height * magnification > 1080 * 2) {
		printf("Error -- combination of image size and magnification too high.\n");
		exit(1);
	}

	if (compare_files_flag) {
		compare_files();
		exit(0);
	}

	create_result_pixbuf();

	create_window_layout(image_width * magnification, image_height * magnification);

	g_idle_add(start_ga_cb, NULL);
	gtk_main();

	exit(0);
}

