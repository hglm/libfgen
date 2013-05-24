/*
    decode.c -- texture decoding functions.

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "decode.h"

static int modifier_table[8][4] = {
	{ 2, 8, -2, -8 },
	{ 5, 17, -5, -17 },
	{ 9, 29, -9, -29 },
	{ 13, 42, -13, -42 },
	{ 18, 60, -18, -60 },
	{ 24, 80, -24, -80 },
	{ 33, 106, -33, -106 },
	{ 47, 183, -47, -183 }
	};

static int clamp(int x) {
	if (x < 0)
		x = 0;
	if (x > 255)
		x = 255;
	return x;
}

static int complement3bitshifted(int x) {
	if (x & 4)
		return ((x & 3) - 4) << 3;	// Note: shift is arithmetic.
	return x << 3;
}

static int complement3bit(int x) {
	if (x & 4)
		return ((x & 3) - 4);
	return x;
}

// Draw a 4x4 pixel block using the ETC1 texture compression data in bitstring. Return 1 is it is a valid block, 0 if
// not (overflow or underflow occurred in differential mode).

int draw_block4x4_etc1(const unsigned char *bitstring, unsigned int *image_buffer, int modes_allowed) {
	int differential_mode = bitstring[3] & 2;
	if (differential_mode) {
		if ((modes_allowed & ETC_MODE_ALLOWED_DIFFERENTIAL) == 0)
			return 0;
	}
	else
		if ((modes_allowed & ETC_MODE_ALLOWED_INDIVIDUAL) == 0)
			return 0;
	int flipbit = bitstring[3] & 1;
	int base_color_subblock1_R;
	int base_color_subblock1_G;
	int base_color_subblock1_B;
	int base_color_subblock2_R;
	int base_color_subblock2_G;
	int base_color_subblock2_B;
	if (differential_mode) {
		base_color_subblock1_R = (bitstring[0] & 0xF8);
		base_color_subblock1_R |= ((base_color_subblock1_R & 224) >> 5);
		base_color_subblock1_G = (bitstring[1] & 0xF8);
		base_color_subblock1_G |= (base_color_subblock1_G & 224) >> 5;
		base_color_subblock1_B = (bitstring[2] & 0xF8);
		base_color_subblock1_B |= (base_color_subblock1_B & 224) >> 5;
		base_color_subblock2_R = (bitstring[0] & 0xF8);				// 5 highest order bits.
		base_color_subblock2_R += complement3bitshifted(bitstring[0] & 7);	// Add difference.
		if (base_color_subblock2_R & 0xFF07)					// Check for overflow.
			return 0;
		base_color_subblock2_R |= (base_color_subblock2_R & 224) >> 5;		// Replicate.
		base_color_subblock2_G = (bitstring[1] & 0xF8);
		base_color_subblock2_G += complement3bitshifted(bitstring[1] & 7);
		if (base_color_subblock2_G & 0xFF07)
			return 0;
		base_color_subblock2_G |= (base_color_subblock2_G & 224) >> 5;
		base_color_subblock2_B = (bitstring[2] & 0xF8);
		base_color_subblock2_B += complement3bitshifted(bitstring[2] & 7);
		if (base_color_subblock2_B & 0xFF07)
			return 0;
		base_color_subblock2_B |= (base_color_subblock2_B & 224) >> 5;
	}
	else {
		base_color_subblock1_R = (bitstring[0] & 0xF0);
		base_color_subblock1_R |= base_color_subblock1_R >> 4;
		base_color_subblock1_G = (bitstring[1] & 0xF0);
		base_color_subblock1_G |= base_color_subblock1_G >> 4;
		base_color_subblock1_B = (bitstring[2] & 0xF0);
		base_color_subblock1_B |= base_color_subblock1_B >> 4;
		base_color_subblock2_R = (bitstring[0] & 0x0F);
		base_color_subblock2_R |= base_color_subblock2_R << 4;
		base_color_subblock2_G = (bitstring[1] & 0x0F);
		base_color_subblock2_G |= base_color_subblock2_G << 4;
		base_color_subblock2_B = (bitstring[2] & 0x0F);
		base_color_subblock2_B |= base_color_subblock2_B << 4;
	}
	int table_codeword1 = (bitstring[3] & 224) >> 5;
	int table_codeword2 = (bitstring[3] & 28) >> 2;
	unsigned int pixel_index_word = ((unsigned int)bitstring[4] << 24) | ((unsigned int)bitstring[5] << 16) |
		((unsigned int)bitstring[6] << 8) | bitstring[7];
	for (int i = 0; i < 16; i++) {
		int pixel_index = ((pixel_index_word & (1 << i)) >> i)		// Least significant bit.
			| ((pixel_index_word & (0x10000 << i)) >> (16 + i - 1));	// Most significant bit.
		int r, g, b;
		if (flipbit == 0) {
			// Two 2x4 blocks side-by-side.
			if (i < 8) {
				// Subblock 1.
				int modifier = modifier_table[table_codeword1][pixel_index];
				r = clamp(base_color_subblock1_R + modifier);
				g = clamp(base_color_subblock1_G + modifier);
				b = clamp(base_color_subblock1_B + modifier);
			}
			else {
				// Subblock 2.
				int modifier = modifier_table[table_codeword2][pixel_index];
				r = clamp(base_color_subblock2_R + modifier);
				g = clamp(base_color_subblock2_G + modifier);
				b = clamp(base_color_subblock2_B + modifier);
			}
		}
		else {
			// Two 4x2 blocks on top of each other.
			if ((i & 2) == 0) {
				// Subblock 1.
				int modifier = modifier_table[table_codeword1][pixel_index];
				r = clamp(base_color_subblock1_R + modifier);
				g = clamp(base_color_subblock1_G + modifier);
				b = clamp(base_color_subblock1_B + modifier);
			}
			else {
				// Subblock 2.
				int modifier = modifier_table[table_codeword2][pixel_index];
				r = clamp(base_color_subblock2_R + modifier);
				g = clamp(base_color_subblock2_G + modifier);
				b = clamp(base_color_subblock2_B + modifier);
			}
		}
		image_buffer[(i & 3) * 4 + ((i & 12) >> 2)] = r + (g << 8) + (b << 16);
	}
	return 1;
}

static int etc2_distance_table[8] = { 3, 6, 11, 16, 23, 32, 41, 64 };

static void draw_block4x4_rgb8_etc2_T_or_H_mode(const unsigned char *bitstring, unsigned int *image_buffer, int mode) {
	int base_color1_R, base_color1_G, base_color1_B;
	int base_color2_R, base_color2_G, base_color2_B;
	int paint_color_R[4], paint_color_G[4], paint_color_B[4];
	int distance;
	if (mode == 0) {
		// T mode.
		base_color1_R = ((bitstring[0] & 0x18) >> 1) | (bitstring[0] & 0x3);
		base_color1_R |= base_color1_R << 4;
		base_color1_G = bitstring[1] & 0xF0;
		base_color1_G |= base_color1_G >> 4;
		base_color1_B = bitstring[1] & 0x0F;
		base_color1_B |= base_color1_B << 4;
		base_color2_R = bitstring[2] & 0xF0;
		base_color2_R |= base_color2_R >> 4;
		base_color2_G = bitstring[2] & 0x0F;
		base_color2_G |= base_color2_G << 4;
		base_color2_B = bitstring[3] & 0xF0;
		base_color2_B |= base_color2_B >> 4;
		// index = (da << 1) | db
		distance = etc2_distance_table[((bitstring[3] & 0x0C) >> 1) | (bitstring[3] & 0x1)];
		paint_color_R[0] = base_color1_R;
		paint_color_G[0] = base_color1_G;
		paint_color_B[0] = base_color1_B;
		paint_color_R[2] = base_color2_R;
		paint_color_G[2] = base_color2_G;
		paint_color_B[2] = base_color2_B;
		paint_color_R[1] = clamp(base_color2_R + distance);
		paint_color_G[1] = clamp(base_color2_G + distance);
		paint_color_B[1] = clamp(base_color2_B + distance);
		paint_color_R[3] = clamp(base_color2_R - distance);
		paint_color_G[3] = clamp(base_color2_G - distance);
		paint_color_B[3] = clamp(base_color2_B - distance);
	}
	else {
		// H mode.
		base_color1_R = (bitstring[0] & 0x78) >> 3;
		base_color1_R |= base_color1_R << 4;
		base_color1_G = ((bitstring[0] & 0x07) << 1) | ((bitstring[1] & 0x10) >> 4);
		base_color1_G |= base_color1_G << 4;
		base_color1_B = (bitstring[1] & 0x08) | ((bitstring[1] & 0x03) << 1) | ((bitstring[2] & 0x80) >> 7);
		base_color1_B |= base_color1_B << 4;
		base_color2_R = (bitstring[2] & 0x78) >> 3;
		base_color2_R |= base_color2_R << 4;
		base_color2_G = ((bitstring[2] & 0x07) << 1) | ((bitstring[3] & 0x80) >> 7);
		base_color2_G |= base_color2_G << 4;
		base_color2_B = (bitstring[3] & 0x78) >> 3;
		base_color2_B |= base_color2_B << 4;
		// da is most significant bit, db is middle bit, least significant bit is
		// (base_color1 value >= base_color2 value).
		int base_color1_value = (base_color1_R << 16) + (base_color1_G << 8) + base_color1_B;
		int base_color2_value = (base_color2_R << 16) + (base_color2_G << 8) + base_color2_B;
		int bit;
		if (base_color1_value >= base_color2_value)
			bit = 1;
		else
			bit = 0;
		distance = etc2_distance_table[(bitstring[3] & 0x04) | ((bitstring[3] & 0x01) << 1) | bit];
		paint_color_R[0] = clamp(base_color1_R + distance);
		paint_color_G[0] = clamp(base_color1_G + distance);
		paint_color_B[0] = clamp(base_color1_B + distance);
		paint_color_R[1] = clamp(base_color1_R - distance);
		paint_color_G[1] = clamp(base_color1_G - distance);
		paint_color_B[1] = clamp(base_color1_B - distance);
		paint_color_R[2] = clamp(base_color2_R + distance);
		paint_color_G[2] = clamp(base_color2_G + distance);
		paint_color_B[2] = clamp(base_color2_B + distance);
		paint_color_R[3] = clamp(base_color2_R - distance);
		paint_color_G[3] = clamp(base_color2_G - distance);
		paint_color_B[3] = clamp(base_color2_B - distance);
	}
	unsigned int pixel_index_word = ((unsigned int)bitstring[4] << 24) | ((unsigned int)bitstring[5] << 16) |
		((unsigned int)bitstring[6] << 8) | bitstring[7];
	for (int i = 0; i < 16; i++) {
		int pixel_index = ((pixel_index_word & (1 << i)) >> i)			// Least significant bit.
			| ((pixel_index_word & (0x10000 << i)) >> (16 + i - 1));	// Most significant bit.
		int r = paint_color_R[pixel_index];
		int g = paint_color_G[pixel_index];
		int b = paint_color_B[pixel_index];
		image_buffer[(i & 3) * 4 + ((i & 12) >> 2)] = r + (g << 8) + (b << 16);
	}

}

static void draw_block4x4_rgb8_etc2_planar_mode(const unsigned char *bitstring, unsigned int *image_buffer) {
	// Each color O, H and V is in 6-7-6 format.
	int RO = (bitstring[0] & 0x7E) >> 1;
	int GO = ((bitstring[0] & 0x1) << 6) | ((bitstring[1] & 0x7E) >> 1);
	int BO = ((bitstring[1] & 0x1) << 5) | (bitstring[2] & 0x18) | ((bitstring[2] & 0x03) << 1) |
		((bitstring[3] & 0x80) >> 7);
	int RH = ((bitstring[3] & 0x7B) >> 1) | (bitstring[3] & 0x1);
	int GH = (bitstring[4] & 0xFE) >> 1;
	int BH = ((bitstring[4] & 0x1) << 5) | ((bitstring[5] & 0xF8) >> 3);
	int RV = ((bitstring[5] & 0x7) << 3) | ((bitstring[6] & 0xE0) >> 5);
	int GV = ((bitstring[6] & 0x1F) << 2) | ((bitstring[7] & 0xC0) >> 6);
	int BV = bitstring[7] & 0x3F;
	RO = (RO << 2) | ((RO & 0x30) >> 4);	// Replicate bits.
	GO = (GO << 1) | ((GO & 0x40) >> 6);
	BO = (BO << 2) | ((BO & 0x30) >> 4);
	RH = (RH << 2) | ((RH & 0x30) >> 4);
	GH = (GH << 1) | ((GH & 0x40) >> 6);
	BH = (BH << 2) | ((BH & 0x30) >> 4);
	RV = (RV << 2) | ((RV & 0x30) >> 4);
	GV = (GV << 1) | ((GV & 0x40) >> 6);
	BV = (BV << 2) | ((BV & 0x30) >> 4);
	for (int y = 0; y < 4; y++)
		for (int x = 0; x < 4; x++) {
			int r = clamp((x * (RH - RO) + y * (RV - RO) + 4 * RO + 2) >> 2);
			int g = clamp((x * (GH - GO) + y * (GV - GO) + 4 * GO + 2) >> 2);
			int b = clamp((x * (BH - BO) + y * (BV - BO) + 4 * BO + 2) >> 2);
			image_buffer[y * 4 + x] = r + (g << 8) + (b << 16);
		}
}

// Draw a 4x4 pixel block using the ETC2 texture compression data in bitstring (format COMPRESSED_RGB8_ETC2).

int draw_block4x4_rgb8_etc2(const unsigned char *bitstring, unsigned int *image_buffer, int modes_allowed) {
	// Figure out the mode.
	if ((bitstring[3] & 2) == 0) {
		// Individual mode.
		return draw_block4x4_etc1(bitstring, image_buffer, modes_allowed);
	}
	if ((modes_allowed & ~ETC_MODE_ALLOWED_INDIVIDUAL) == 0)
		return 0;
#if 0
	int R = (bitstring[0] & 0xF8) >> 3;;
	R += complement3bit(bitstring[0] & 7);
	int G = (bitstring[1] & 0xF8) >> 3;
	G += complement3bit(bitstring[1] & 7);
	int B = (bitstring[2] & 0xF8) >> 3;
	B += complement3bit(bitstring[2] & 7);
	if (R < 0 || R > 31) {
		// T mode.
		if ((modes_allowed & ETC2_MODE_ALLOWED_T) == 0)
			return 0;
		draw_block4x4_rgb8_etc2_T_or_H_mode(bitstring, image_buffer, 0);
		return 1;
	}
	else
	if (G < 0 || G > 31) {
		// H mode.
		if ((modes_allowed & ETC2_MODE_ALLOWED_H) == 0)
			return 0;
		draw_block4x4_rgb8_etc2_T_or_H_mode(bitstring, image_buffer, 1);
		return 1;
	}
	else
	if (B < 0 || B > 31) {
		// Planar mode.
		if ((modes_allowed & ETC2_MODE_ALLOWED_PLANAR) == 0)
			return 0;
		draw_block4x4_rgb8_etc2_planar_mode(bitstring, image_buffer);
		return 1;
	}
	else {
		// Differential mode.
		return draw_block4x4_etc1(bitstring, image_buffer, modes_allowed);
	}

#else
	int R = (bitstring[0] & 0xF8);
	R += complement3bitshifted(bitstring[0] & 7);
	int G = (bitstring[1] & 0xF8);
	G += complement3bitshifted(bitstring[1] & 7);
	int B = (bitstring[2] & 0xF8);
	B += complement3bitshifted(bitstring[2] & 7);
	if (R & 0xFF07) {
		// T mode.
		if ((modes_allowed & ETC2_MODE_ALLOWED_T) == 0)
			return 0;
		draw_block4x4_rgb8_etc2_T_or_H_mode(bitstring, image_buffer, 0);
		return 1;
	}
	else
	if (G & 0xFF07) {
		// H mode.
		if ((modes_allowed & ETC2_MODE_ALLOWED_H) == 0)
			return 0;
		draw_block4x4_rgb8_etc2_T_or_H_mode(bitstring, image_buffer, 1);
		return 1;
	}
	else
	if (B & 0xFF07) {
		// Planar mode.
		if ((modes_allowed & ETC2_MODE_ALLOWED_PLANAR) == 0)
			return 0;
		draw_block4x4_rgb8_etc2_planar_mode(bitstring, image_buffer);
		return 1;
	}
	else {
		// Differential mode.
		return draw_block4x4_etc1(bitstring, image_buffer, modes_allowed);
	}
#endif
}

// Only determine the mode of the ETC2 RGB8 block.
// Return one of the following values:
// 0	Individual mode
// 1	Differential mode
// 2	T mode
// 3	H mode
// 4	Planar mode	

int block4x4_rgb8_etc2_get_mode(const unsigned char *bitstring) {
	// Figure out the mode.
	if ((bitstring[3] & 2) == 0)
		// Individual mode.
		return 0;
#if 0
	int R = (bitstring[0] & 0xF8) >> 3;;
	R += complement3bit(bitstring[0] & 7);
	int G = (bitstring[1] & 0xF8) >> 3;
	G += complement3bit(bitstring[1] & 7);
	int B = (bitstring[2] & 0xF8) >> 3;
	B += complement3bit(bitstring[2] & 7);
	if (R < 0 || R > 31)
		// T mode.
		return 2;
	else
	if (G < 0 || G > 31)
		// H mode.
		return 3;
	else
	if (B < 0 || B > 31)
		// Planar mode.
		return 4;
	else
		// Differential mode.
		return 1;
#else
	int R = (bitstring[0] & 0xF8);
	R += complement3bitshifted(bitstring[0] & 7);
	int G = (bitstring[1] & 0xF8);
	G += complement3bitshifted(bitstring[1] & 7);
	int B = (bitstring[2] & 0xF8);
	B += complement3bitshifted(bitstring[2] & 7);
	if (R & 0xFF07)
		// T mode.
		return 2;
	else
	if (G & 0xFF07)
		// H mode.
		return 3;
	else
	if (B & 0xFF07)
		// Planar mode.
		return 4;
	else
		// Differential mode.
		return 1;
#endif
}

// Draw a 4x4 pixel block using the DXT1 texture compression data in bitstring.

int draw_block4x4_dxt1(const unsigned char *bitstring, unsigned int *image_buffer) {
	unsigned int colors = *(unsigned int *)&bitstring[0];
	// Decode the two 5-6-5 RGB colors.
	int color_r[4], color_g[4], color_b[4];
	color_b[0] = (colors & 0x0000001F) << 3;
	color_g[0] = (colors & 0x000007E0) >> (5 - 2);
	color_r[0] = (colors & 0x0000F800) >> (11 - 3);
	color_b[1] = (colors & 0x001F0000) >> (16 - 3);
	color_g[1] = (colors & 0x07E00000) >> (21 - 2);
	color_r[1] = (colors & 0xF8000000) >> (27 - 3);
	if ((colors & 0xFFFF) > ((colors & 0xFFFF0000) >> 16)) {
		color_r[2] = (2 * color_r[0] + color_r[1]) / 3;
		color_g[2] = (2 * color_g[0] + color_g[1]) / 3;
		color_b[2] = (2 * color_b[0] + color_b[1]) / 3;
		color_r[3] = (color_r[0] + 2 * color_r[1]) / 3;
		color_g[3] = (color_g[0] + 2 * color_g[1]) / 3;
		color_b[3] = (color_b[0] + 2 * color_b[1]) / 3;
	}
	else {
		color_r[2] = (color_r[0] + color_r[1]) / 2;
		color_g[2] = (color_g[0] + color_g[1]) / 2;
		color_b[2] = (color_b[0] + color_b[1]) / 2;
		color_r[3] = color_g[3] = color_b[3] = 0;
	}
	unsigned int pixels = *(unsigned int *)&bitstring[4];
	for (int i = 0; i < 16; i++) {
		int pixel = (pixels >> (i * 2)) & 0x3;
		image_buffer[i] = color_r[pixel] + (color_g[pixel] << 8) + (color_b[pixel] << 16);
	}
	return 1;
}


