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

// bil_block.c
// Wes Kendall
// 07/18/2010
// Functions for processing blocks in BIL. These functions include the
// extraction of an N-d block from a larger block and functions for the
// comparison of blocks when sorting. Some other functions are defined that
// add a block to the internal structure of BIL.
////

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "bil.h"
#include "bil_block.h"
#include "bil_misc.h"

// BIL_Block_add
// Adds a block to BIL's internal structure and returns the pointer to it.
////
BIL_Block* BIL_Block_add() {
  BIL->num_blocks++;
  BIL->blocks =
    (BIL_Block*)BIL_Misc_realloc(BIL->blocks,
                                 sizeof(BIL_Block) * BIL->num_blocks);
  return &(BIL->blocks[BIL->num_blocks - 1]);
}

// BIL_Block_init
// Initializes values of a block and assigns the block an id.
////
void BIL_Block_init(BIL_Block* block, int num_dims, const int* block_starts,
                    const int* block_sizes, const char* file_name,
                    void** buffer) {
  memset(block, 0, sizeof(BIL_Block));
  block->num_dims = num_dims;
  memcpy(block->starts, block_starts, num_dims * sizeof(int));
  memcpy(block->sizes, block_sizes, num_dims * sizeof(int));
  strcpy(block->file_name, file_name);
  block->request_rank = BIL->world_rank;
  block->io_group = -1;

  // Find the total size of the block.
  int i;
  block->total_size = 1;
  for (i = 0; i < block->num_dims; i++)  {
    block->total_size *= block->sizes[i];
  }

  // Assign a unique id to the block.
  block->block_id = BIL->num_blocks - 1;

  // Assign the buffer where data where read/write data is stored.
  block->data = *buffer;
  block->user_data_buffer = buffer;
}

// BIL_Block_extract_4d
// Extracts a 4d subset from a block.
////
int BIL_Block_extract_4d(const void *data, const int *data_start,
                         const int *data_size, void *block,
                         const int *block_start, const int *block_size,
                         int var_size) {   
  int data_offset = 0, i, j, k;
  for (i = block_start[0]; i < block_start[0] + block_size[0]; i++) {
    for (j = block_start[1]; j < block_start[1] + block_size[1]; j++) {
      for (k = block_start[2]; k < block_start[2] + block_size[2]; k++) {
        memcpy(block + data_offset, 
               data + ((((i - data_start[0]) * data_size[1] 
                 * data_size[2] * data_size[3])
                 + ((j - data_start[1]) * data_size[2] * data_size[3])
                 + ((k - data_start[2]) * data_size[3])
                 + (block_start[3] - data_start[3])) * var_size),
               block_size[3]  * var_size);
        data_offset += block_size[3] * var_size;
      }
    }
  }
  return data_offset;
}

// BIL_Block_extract_3d
// Extracts a 3d subset from a block.
////
int BIL_Block_extract_3d(const void *data, const int *data_start,
                         const int *data_size, void *block,
                         const int *block_start, const int *block_size,
                         int var_size) {   
  int data_offset = 0, i, j;
  for (i = block_start[0]; i < block_start[0] + block_size[0]; i++) {
    for (j = block_start[1]; j < block_start[1] + block_size[1]; j++) {
      memcpy(block + data_offset,
             data + ((((i - data_start[0]) * data_size[1] * data_size[2])
               + ((j - data_start[1]) * data_size[2])
               + (block_start[2] - data_start[2])) * var_size),
             block_size[2]  * var_size);
      data_offset += block_size[2] * var_size;
    }
  }
  return data_offset;
}

// BIL_Block_extract_2d
// Extracts a 2d subset from a block.
////
int BIL_Block_extract_2d(const void *data, const int *data_start,
                         const int *data_size, void *block,
                         const int *block_start, const int *block_size,
                         int var_size) {   
  int data_offset = 0, i;
  for (i = block_start[0]; i < block_start[0] + block_size[0]; i++) {
    memcpy(block + data_offset,
           data + ((((i - data_start[0]) * data_size[1])
             + (block_start[1] - data_start[1])) * var_size),
           block_size[1]  * var_size);
    data_offset += block_size[1] * var_size;
  }
  return data_offset;
}

// BIL_Block_extract_1d
// Extracts a 1d subset from a block.
////
int BIL_Block_extract_1d(const void *data, const int *data_start,
                         const int *data_size, void *block,
                         const int *block_start, const int *block_size,
                         int var_size) {   
  memcpy(block, data + ((block_start[0] - data_start[0]) * var_size),
         block_size[0] * var_size);
  return block_size[0] * var_size;
}
