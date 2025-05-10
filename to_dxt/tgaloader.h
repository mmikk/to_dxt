/**
 *  Copyright (C) 2011 by Morten S. Mikkelsen
 *
 *  This software is provided 'as-is', without any express or implied
 *  warranty.  In no event will the authors be held liable for any damages
 *  arising from the use of this software.
 *
 *  Permission is granted to anyone to use this software for any purpose,
 *  including commercial applications, and to alter it and redistribute it
 *  freely, subject to the following restrictions:
 *
 *  1. The origin of this software must not be misrepresented; you must not
 *     claim that you wrote the original software. If you use this software
 *     in a product, an acknowledgment in the product documentation would be
 *     appreciated but is not required.
 *  2. Altered source versions must be plainly marked as such, and must not be
 *     misrepresented as being the original software.
 *  3. This notice may not be removed or altered from any source distribution.
 */

#ifndef __TGALOADER_H__
#define __TGALOADER_H__

#include <stdio.h>

// bitmaps are received and given in the traditional zig zag pattern.
// first 3 bytes represent the upper-left pixel and the last 3 bytes
// represent the lower-right pixel. The pixels are ordered as r, g, b, r, ...

bool ReadTGA(unsigned char ** ppbitmap, int * piWidth, int * piHeight, int * piNrChannels, FILE * fptr_in);
void WriteTGA_24Bit(FILE * fptr_out, const unsigned char bitmap[], const int iWidth, const int iHeight);
void WriteTGA_32Bit(FILE * fptr_out, const unsigned char bitmap[], const int iWidth, const int iHeight);


#endif