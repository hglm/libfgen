/*
    decode.h -- part of texgenpack, a texture compressor using fgen.

    fgen -- Library for optimization using a genetic algorithm or particle swarm optimization.
    Copyright 2013, Harm Hanemaaijer

    This file is part of texgenpack.

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

#define ETC_MODE_ALLOWED_INDIVIDUAL	1
#define ETC_MODE_ALLOWED_DIFFERENTIAL	2
#define ETC2_MODE_ALLOWED_T		4
#define ETC2_MODE_ALLOWED_H		8
#define ETC2_MODE_ALLOWED_PLANAR	16
#define ETC_MODE_ALLOWED_ALL		3
#define ETC2_MODE_ALLOWED_ALL		31

// Draw (decompress) a 64-bit 4x4 pixel block.
int draw_block4x4_etc1(const unsigned char *bitstring, unsigned int *image_buffer, int modes_allowed);
int draw_block4x4_rgb8_etc2(const unsigned char *bitstring, unsigned int *image_buffer, int modes_allowed);
int draw_block4x4_dxt1(const unsigned char *bitstring, unsigned int *image_buffer);

// Return ETC2 mode number from 0 to 4.
int block4x4_rgb8_etc2_get_mode(const unsigned char *bitstring);

