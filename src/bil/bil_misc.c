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

// bil_misc.c
// Wes Kendall
// 07/18/2010
// Miscellaneous functions for BIL, including memory allocation and prefix sum
// computation.
////

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <errno.h>

#include "bil.h"
#include "bil_misc.h"

// BIL_Misc_error
// Prints an error message and aborts.
////
void BIL_Misc_error(BIL_Error_Type error_code) {
  fprintf(stderr, "BIL internal error %d\n", error_code);
  MPI_Abort(BIL->world_comm, 1);
}

void* BIL_Misc_malloc(size_t size) {
  if (size == 0) {
    return NULL;
  }
  void* buf = malloc(size);
  if (buf == NULL) {
    fprintf(stderr, "BIL malloc error for size %d, %lf MB\n",
            (int)size, size / (1024.0 * 1024.0));
  }
  assert(buf != NULL);
  return buf;
}

// BIL_Misc_malloc_z
// A wrapper around malloc that checks for NULL return value and also zeroes out
// the memory.
////
void* BIL_Misc_malloc_z(size_t size) {
  if (size == 0) {
    return NULL;
  }
  void* buf = malloc(size);
  assert(buf != NULL);
  memset(buf, 0, size);
  return buf;
}

// BIL_Misc_realloc
// A wrapper around realloc that checks for NULL return.
////
void* BIL_Misc_realloc(void *orig_buf, size_t size) {
  void* buf = realloc(orig_buf, size);
  assert(buf != NULL);
  return buf;
}

// BIL_Misc_free
// A wrapper around free that checks for NULL values.
////
void BIL_Misc_free(void *buf) {
  if (buf != NULL) {
    free(buf);
  }
}

// BIL_Misc_prefix_sum
// Given an array of sizes, produce an array of offsets.
////
void BIL_Misc_prefix_sum(int num_items, const int* sizes, int* offsets) {
  int i;
  offsets[0] = 0;
  for (i = 1; i < num_items; i++) {
    offsets[i] = sizes[i - 1] + offsets[i - 1];
  }
}

// BIL_Misc_array_sum
// Given an array of integers, return the sum.
////
int BIL_Misc_array_sum(const int* array, int num_elements) {
  int i, sum;
  for (sum = 0, i = 0; i < num_elements; i++) {
    sum += array[i];
  }
  return sum;
}

// BIL_Misc_gatherv
// Wrapper around MPI_Gatherv that computes displacements from the given counts.
////
void BIL_Misc_gatherv(void* send_data, int send_count, MPI_Datatype send_type,
                      void* recv_data, int* recv_counts,
                      MPI_Datatype recv_type) {
  int* recv_displs = NULL;
  if (BIL->world_rank == 0) { 
    recv_displs = BIL_Misc_malloc(sizeof(int) * BIL->world_size);
    BIL_Misc_prefix_sum(BIL->world_size, recv_counts, recv_displs);
  }
  MPI_Gatherv(send_data, send_count, send_type, recv_data, recv_counts,
              recv_displs, recv_type, 0, BIL->world_comm);
  if (BIL->world_rank == 0) { 
    BIL_Misc_free(recv_displs);
  }
}

// BIL_Misc_scatterv
// Wrapper around MPI_Scatterv that computes displacements from the given
// counts.
////
void BIL_Misc_scatterv(void* send_data, int* send_counts,
                       MPI_Datatype send_type, void* recv_data, int recv_count,
                       MPI_Datatype recv_type) {
  int* send_displs = NULL;
  if (BIL->world_rank == 0) { 
    send_displs = BIL_Misc_malloc(sizeof(int) * BIL->world_size);
    BIL_Misc_prefix_sum(BIL->world_size, send_counts, send_displs);
  }
  MPI_Scatterv(send_data, send_counts, send_displs, send_type, recv_data,
               recv_count, recv_type, 0, BIL->world_comm);
  if (BIL->world_rank == 0) { 
    BIL_Misc_free(send_displs);
  }
}
