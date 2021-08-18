/***************************************************************************
 *   Copyright (C) 2010  Wes Kendall                                       *
 *   kendall@eecs.utk.edu                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 3 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU Lesser General Public License for more details.                   *
 *                                                                         *
 *   You should have received a copy of the GNU Lesser General Public      *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

// bil_block.h
// Wes Kendall
// 07/18/2010
// Prototypes for block processing functions in BIL.
////

#ifndef __BIL_BLOCK_H
#define __BIL_BLOCK_H

#include "bil.h"

#if HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

BIL_Block* BIL_Block_add();

void BIL_Block_init(BIL_Block* block, int num_dims, const int* block_start,
                    const int* block_size, const char* file_name,
                    void** buffer);

int BIL_Block_extract_4d(const void *data, const int *data_start,
                         const int *data_size, void *block,
                         const int *block_start, const int *block_size,
                         int var_size);

int BIL_Block_extract_3d(const void *data, const int *data_start,
                         const int *data_size, void *block,
                         const int *block_start, const int *block_size,
                         int var_size);

int BIL_Block_extract_2d(const void *data, const int *data_start,
                         const int *data_size, void *block,
                         const int *block_start, const int *block_size,
                         int var_size);

int BIL_Block_extract_1d(const void *data, const int *data_start,
                         const int *data_size, void *block,
                         const int *block_start, const int *block_size,
                         int var_size);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // __BIL_BLOCK_H
