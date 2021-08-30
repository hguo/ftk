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

// bil_misc.h
// Wes Kendall
// 07/18/2010
// The function prototypes for the utility functions in BIL.
////

#ifndef __BIL_MISC_H
#define __BIL_MISC_H

#include <errno.h>
#include <stdio.h>

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

void BIL_Misc_error(BIL_Error_Type error_code);

void* BIL_Misc_malloc(size_t size);

void* BIL_Misc_malloc_z(size_t size);

void* BIL_Misc_realloc(void* orig_buf, size_t size);

void BIL_Misc_free(void* buf);

void BIL_Misc_prefix_sum(int num_items, const int* sizes, int* offsets);

int BIL_Misc_array_sum(const int* array, int num_elements);

void BIL_Misc_gatherv(void* send_data, int send_count, MPI_Datatype send_type,
                      void* recv_data, int* recv_counts,
                      MPI_Datatype recv_type);

void BIL_Misc_scatterv(void* send_data, int* send_counts,
                       MPI_Datatype send_type, void* recv_data, int recv_count,
                       MPI_Datatype recv_type);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // __BIL_MISC_H
