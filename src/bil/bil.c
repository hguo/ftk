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

// bil.c
// Wes Kendall
// 07/18/2010
// The functions for the main BIL interface and the only ones available to the
// user.
////

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>
#ifndef DISABLE_PNETCDF
#include <pnetcdf.h>
#endif // DISABLE_PNETCDF

#include "bil.h"
#include "bil_misc.h"
#include "bil_block.h"
#include "bil_timing.h"
#include "bil_sched.h"
#include "bil_pio.h"

// checks for initialization of BIL
#define BIL_INIT_CHECK assert(BIL != NULL);

// the main global structure of BIL
BIL_Info* BIL = NULL;

void BIL_Read() {
  BIL_INIT_CHECK;

  BIL_Timing_total_start();

  int* blocks_to_group_map;
  int num_groups;
  BIL_Sched_IO io_sched, inv_io_sched;
  BIL_Sched_get(&io_sched, &inv_io_sched, &num_groups, &blocks_to_group_map);

  BIL_Pio_issue(&io_sched, &inv_io_sched, num_groups, blocks_to_group_map);

  // The data is now stored in BIL's global struct. The user provided void**
  // buffers that either had a buffer already allocated or had a NULL buffer
  // that BIL allocated for them. For the case where BIL allocated the buffer,
  // go through the blocks and set the user provided buffer to the allocated
  // one.
#if 0 // commented out by hanqi
  int i;
  for (i = 0; i < BIL->num_blocks; i++) {
    *(BIL->blocks[i].user_data_buffer) = BIL->blocks[i].data;
  }
#endif
  BIL_Timing_total_stop();
 
  BIL_Misc_free(io_sched.io_blocks);
  BIL_Misc_free(inv_io_sched.io_blocks);
  BIL_Misc_free(blocks_to_group_map);
  BIL_Misc_free(BIL->blocks);
  BIL->num_blocks = 0;
  BIL->blocks = NULL;
}

void BIL_Add_block_raw(int num_dims, const int *data_dims,
                       const int *block_start, const int *block_size,
                       const char *file_name, MPI_Datatype var_type,
                       void** buffer) {
  BIL_INIT_CHECK;
  assert(strlen(file_name) < BIL_MAX_FILE_NAME_SIZE);

  // Add a block to BIL's internal structure.
  BIL_Block *block = BIL_Block_add();

  // Initialize values general to every block.
  BIL_Block_init(block, num_dims, block_start, block_size, file_name,
                 buffer);

  // Initialize values specific to raw blocks.
  memcpy(block->file_dim_sizes, data_dims, num_dims * sizeof(int));  
  MPI_Type_size(var_type, &(block->var_size));
  block->io_type = BIL_RAW;

  return;
}

#ifndef DISABLE_PNETCDF
void BIL_Add_block_nc(int num_dims, const int *block_start,
                      const int *block_size, const char *file_name,
                      const char *var_name, void** buffer) {
  BIL_INIT_CHECK;
  assert(strlen(file_name) < BIL_MAX_FILE_NAME_SIZE);
  assert(strlen(var_name) < BIL_MAX_VAR_NAME_SIZE);
  
  // Add a block to BIL's internal structure.
  BIL_Block* block = BIL_Block_add();

  // Initialize values general to every block.
  BIL_Block_init(block, num_dims, block_start, block_size, file_name, buffer);

  // Initialize values specific to netcdf I/O.
  strcpy(block->var_name, var_name);
  block->io_type = BIL_NC;
}
#endif // DISABLE_PNETCDF

void BIL_Set_io_hints(MPI_Info io_hints) {
  // Free any io hints that might already be set.
  if (BIL->io_hints != MPI_INFO_NULL) {
    MPI_Info_free(&(BIL->io_hints));
  }
  if (io_hints != MPI_INFO_NULL) {
    assert(MPI_Info_dup(io_hints, &(BIL->io_hints)) == MPI_SUCCESS);
  }
}

void BIL_Set_io_header_size(int io_header_size) {
  BIL->io_header_size = io_header_size;
}

void BIL_Init(const MPI_Comm world_comm) {
  assert(BIL == NULL);

  BIL = (BIL_Info*)BIL_Misc_malloc(sizeof(BIL_Info));
  BIL->world_comm = world_comm;
  MPI_Comm_rank(BIL->world_comm, &(BIL->world_rank));
  MPI_Comm_size(BIL->world_comm, &(BIL->world_size));
  BIL->num_blocks = 0;
  BIL->blocks = NULL;

  MPI_Type_contiguous(sizeof(BIL_Block), MPI_BYTE, &(BIL->bil_block_type));
  MPI_Type_commit(&(BIL->bil_block_type));

  BIL_Timing_init(&BIL->timing);

  // Set the user defined parameters to their NULL values.
  BIL->io_hints = MPI_INFO_NULL;
  BIL->io_header_size = 0;
}

void BIL_Finalize() {
  BIL_INIT_CHECK;

  if (BIL->num_blocks > 0) {
    BIL_Misc_free(BIL->blocks);
  }
  
  MPI_Type_free(&(BIL->bil_block_type));
  BIL_Misc_free(BIL);
  BIL = NULL;
}
