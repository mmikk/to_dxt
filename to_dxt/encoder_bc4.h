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

#ifndef __ENCODERBC4_H__
#define __ENCODERBC4_H__

#include "miniutils.h"

// block_in is 16 entries long of uint16 single channel
// uResultingPixels[] is the visual block result after reduction
// block_out is the actual bc4 encoding
void ConvertBlockToBC4(const uint16 block_in[], uint64 block_out[], uint16 uResultingPixels[]);
void DecodeBC4(uint16 uResultingPixels[], const uint64 block_bc4);


#endif
