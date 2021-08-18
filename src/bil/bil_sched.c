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

// bil_sched.c
// Wes Kendall
// 07/18/2010
// Functions for scheduling blocks to processes in BIL.
////

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>

#include "bil.h"
#include "bil_sched.h"
#include "bil_timing.h"
#include "bil_misc.h"
#include "bil_block.h"

int BIL_Sched_compare_file_and_var_name(const void* a, const void* b) {
  BIL_Block* block_a = *(BIL_Block**)a;
  BIL_Block* block_b = *(BIL_Block**)b;
  int file_name_cmp = memcmp(block_a->file_name, block_b->file_name,
                             BIL_MAX_FILE_NAME_SIZE);
  if (file_name_cmp < 0) {
    return -1;
  } else if (file_name_cmp > 0) {
    return 1;
  } else {
    return memcmp(block_a->var_name, block_b->var_name, BIL_MAX_VAR_NAME_SIZE);
  }
}

int BIL_Sched_compare_group_name(const void* a, const void* b) {
  return memcmp(a, b, BIL_MAX_FILE_NAME_SIZE + BIL_MAX_VAR_NAME_SIZE);
}

void BIL_Sched_merge_unique_group_names(char** group_names_ret, 
                                        int* num_unique_groups_ret) {
  char* group_names = *group_names_ret;
  int num_groups = *num_unique_groups_ret;

  int group_name_size = BIL_MAX_FILE_NAME_SIZE + BIL_MAX_VAR_NAME_SIZE;
  qsort(group_names, num_groups, group_name_size, BIL_Sched_compare_group_name);

  int i;
  int num_unique_groups = 1;
  for (i = 1; i < num_groups; i++) {
    if (memcmp(group_names + (i * group_name_size),
               group_names + ((i - 1) * group_name_size), 
               group_name_size) != 0) {
      memmove(group_names + (num_unique_groups * group_name_size),
              group_names + (i * group_name_size), group_name_size);
      num_unique_groups++;
    }
  }
  *group_names_ret = group_names;
  *num_unique_groups_ret = num_unique_groups;
  return;
}

void BIL_Sched_swap_and_merge_unique_group_names(
    int is_sender, int swapper_rank, char** unique_group_names_ret,
    int* num_unique_groups_ret) {
  char* unique_group_names = *unique_group_names_ret;
  int num_unique_groups = *num_unique_groups_ret;

  int group_name_size = BIL_MAX_FILE_NAME_SIZE + BIL_MAX_VAR_NAME_SIZE;
  if (is_sender == 1) {
    MPI_Send((char*)unique_group_names, group_name_size * num_unique_groups,
             MPI_CHAR, swapper_rank, 0, BIL->world_comm);
  } else {
    // Probe for the message size.
    MPI_Status status;
    MPI_Probe(swapper_rank, 0, BIL->world_comm, &status);
    int group_name_recv_amount;
    MPI_Get_count(&status, MPI_CHAR, &group_name_recv_amount);
    unique_group_names = BIL_Misc_realloc(unique_group_names,
                                          num_unique_groups * group_name_size +
                                          group_name_recv_amount);
    // Receive group names from the outer processes.
    MPI_Recv(unique_group_names + (num_unique_groups * group_name_size),
             group_name_recv_amount, MPI_CHAR,
             swapper_rank, 0, BIL->world_comm, MPI_STATUS_IGNORE);

    num_unique_groups += group_name_recv_amount / group_name_size;
    // Merge the group names into an array of unique names.
    BIL_Sched_merge_unique_group_names(&unique_group_names, 
                                       &num_unique_groups);
    *num_unique_groups_ret = num_unique_groups;
    *unique_group_names_ret = unique_group_names;
  }
}

void BIL_Sched_gather_group_names(char** unique_group_names_ret,
                                  int* num_unique_groups_ret) {
  BIL_Block** sorted_groups =
    BIL_Misc_malloc(sizeof(BIL_Block*) * BIL->num_blocks);

  int i;
  for (i = 0; i < BIL->num_blocks; i++) {
    sorted_groups[i] = &(BIL->blocks[i]);
  }
  qsort(sorted_groups, BIL->num_blocks,
        sizeof(BIL_Block*), BIL_Sched_compare_file_and_var_name);
  
  // Find out how many unique groups you have.
  int num_unique_groups = 1;
  for (i = 1; i < BIL->num_blocks; i++) {
    if (memcmp(sorted_groups[i]->file_name, sorted_groups[i - 1]->file_name,
               BIL_MAX_FILE_NAME_SIZE) != 0 ||
        memcmp(sorted_groups[i]->var_name, sorted_groups[i - 1]->var_name,
               BIL_MAX_VAR_NAME_SIZE) != 0) {
      // Move the unique group to the beginning of the sorted groups list.
      memmove(sorted_groups + num_unique_groups, sorted_groups + i,
              sizeof(BIL_Block*));
      num_unique_groups++;
    }
  }

  // Now the sorted_groups contains all of the unique group names in order.
  // Lets send these babies to each other with some MPI madness. This is
  // perfomed by an algorithm similar to merge sort, with the exception that
  // the sorting only keeps unique values.
  // First, copy all your unique group names into a single array.
  int group_name_size = BIL_MAX_FILE_NAME_SIZE + BIL_MAX_VAR_NAME_SIZE;
  char* unique_group_names = 
    BIL_Misc_malloc(group_name_size * num_unique_groups);
  for (i = 0; i < num_unique_groups; i++) {
    memcpy(unique_group_names + (group_name_size * i),
           sorted_groups[i]->file_name, BIL_MAX_FILE_NAME_SIZE);
    memcpy(unique_group_names + (group_name_size * i) + BIL_MAX_FILE_NAME_SIZE,
           sorted_groups[i]->var_name, BIL_MAX_VAR_NAME_SIZE);
  }
  BIL_Misc_free(sorted_groups);

  // Find the lowest power of two process count.
  int lowest_power_of_two = 1;
  while (lowest_power_of_two < BIL->world_size) {
    lowest_power_of_two <<= 1;
  }

  // Handle the case where there is a non-power-of-two process count.
  // Send all of the unique group names from the outer processes to the
  // processes that are inside the first group that is a power of two.
  if (lowest_power_of_two != BIL->world_size) {
    lowest_power_of_two >>= 1;
    if (BIL->world_rank >= lowest_power_of_two) {
      // Send your group names to the inner processes.
      BIL_Sched_swap_and_merge_unique_group_names(1,
          BIL->world_rank - lowest_power_of_two, &unique_group_names,
          &num_unique_groups);
    } else if (BIL->world_rank < BIL->world_size - lowest_power_of_two) {
      BIL_Sched_swap_and_merge_unique_group_names(0,
          BIL->world_rank + lowest_power_of_two, &unique_group_names,
          &num_unique_groups);
    }
  }
  // Now we have all the data held by a non-power-of-two process count. Do the
  // merging.
  while (BIL->world_rank < lowest_power_of_two && lowest_power_of_two != 1) {
    if (BIL->world_rank >= lowest_power_of_two / 2) {
      BIL_Sched_swap_and_merge_unique_group_names(1,
          BIL->world_rank - (lowest_power_of_two / 2), &unique_group_names,
          &num_unique_groups);
    } else {
      BIL_Sched_swap_and_merge_unique_group_names(0,
          BIL->world_rank + (lowest_power_of_two / 2), &unique_group_names,
          &num_unique_groups);
    }
    lowest_power_of_two >>= 1;
  }
  MPI_Barrier(BIL->world_comm);

  // Now process 0 has all the unique group names. Broadcast the names to
  // everyone.
  MPI_Bcast(&num_unique_groups, 1, MPI_INT, 0, BIL->world_comm);
  if (BIL->world_rank != 0) {
    unique_group_names = BIL_Misc_realloc(unique_group_names,
                                          num_unique_groups * group_name_size);
  }
  MPI_Bcast(unique_group_names, num_unique_groups * group_name_size, MPI_CHAR,
            0, BIL->world_comm);

  *unique_group_names_ret = unique_group_names;
  *num_unique_groups_ret = num_unique_groups;
}
    
void BIL_Sched_get_io_block_assignment(int num_dims,
                                       int io_group, int io_group_rank,
                                       int num_procs_per_group, int* bb_min,
                                       int* bb_max, int* file_res, 
                                       char* group_name, int var_size,
                                       BIL_Block* io_block_ret) {
  // If the io group rank is -1, compute the block io assignment for all ranks
  // within the group.
  int compute_all_block_assignment = (io_group_rank == -1) ? 1 : 0;
  if (compute_all_block_assignment) {
    // To compute all the block assignments, set the group rank to the last
    // one.
    io_group_rank = num_procs_per_group - 1;
  }

  // Compute the total size of the request using the bounding box information.
  int io_total_size = 1;
  int io_size[num_dims];
  int j, k;
  for (k = 0; k < num_dims; k++) {
    io_size[k] = bb_max[k] - bb_min[k];
    io_total_size *= io_size[k];
  }

  int my_io_size[num_dims];
  int io_assigned[num_dims];
  memset(io_assigned, 0, sizeof(int) * num_dims);

  int which_dim_to_start_assigning = 0;
  int num_procs_to_assign = 0;
  for (j = 0; j <= io_group_rank; j++) {
    if (num_procs_to_assign == 0) {
      for (k = num_dims - 1; k >= 0; k--) {
        if (io_assigned[k] != 0) {
          break;
        }
      }
      // Initial case when no I/O has been assigned.
      k = (k == -1) ? 0 : k;

      which_dim_to_start_assigning = k;
      num_procs_to_assign = num_procs_per_group - j;
      for (k = 0; k < which_dim_to_start_assigning; k++) {
        num_procs_to_assign /= (io_size[k] - io_assigned[k]);
      }
      assert(num_procs_to_assign != 0);
 
      for (k = which_dim_to_start_assigning; k < num_dims; k++) {
        if (num_procs_to_assign > io_size[k] - io_assigned[k]) {
          num_procs_to_assign /= io_size[k] - io_assigned[k];
        } else {
          break;
        }
      }
      assert(k < num_dims);
    }

    memset(my_io_size, 0, sizeof(int) * num_dims);
    my_io_size[k] = (io_size[k] - io_assigned[k]) / num_procs_to_assign;
    // Adjust the my_io_size array appropriately.
    int l;
    for (l = 0; l < k; l++) {
      my_io_size[l] = 1;
    }
    for (l = k + 1; l < num_dims; l++) {
      my_io_size[l] = io_size[l];
    }
    if (compute_all_block_assignment == 1 || j == io_group_rank) {
      int which_block = (compute_all_block_assignment == 1) ? j : 0;
      // Copy return value(s)
      io_block_ret[which_block].num_dims = num_dims;
      io_block_ret[which_block].io_group = io_group;
      io_block_ret[which_block].var_size = var_size;
      memcpy(io_block_ret[which_block].file_dim_sizes, file_res, 
             sizeof(int) * BIL_MAX_NUM_DIMS);
      memcpy(io_block_ret[which_block].starts, io_assigned,
             sizeof(int) * BIL_MAX_NUM_DIMS);
      // Offset the starts by the bounding box mins
      for (l = 0; l < num_dims; l++) {
        io_block_ret[which_block].starts[l] += bb_min[l];
      }
      memcpy(io_block_ret[which_block].sizes, my_io_size,
             sizeof(int) * BIL_MAX_NUM_DIMS);
      memcpy(io_block_ret[which_block].file_name, group_name,
             BIL_MAX_FILE_NAME_SIZE);
      memcpy(io_block_ret[which_block].var_name,
             group_name + BIL_MAX_FILE_NAME_SIZE, BIL_MAX_VAR_NAME_SIZE);
    }
    // Adjust the io_assigned array appropriately
    io_assigned[k] += my_io_size[k];

    num_procs_to_assign--;
    while (io_assigned[k] == io_size[k]) {
      io_assigned[k] = 0;
      k--;
      if (k >= 0) {
        io_assigned[k]++;
      } else {
        break;
      }
    }
  }

  if (compute_all_block_assignment == 1) {
    for (k = 0; k < num_dims; k++) {
      assert(io_assigned[k] == 0);
    }
  }
}

int BIL_Sched_intersect_block(const BIL_Block* block_a,
                              const BIL_Block* block_b,
                              BIL_Block* ret_block) {
  int i;
  BIL_Block intersect_block;
  for (i = 0; i < block_a->num_dims; i++) {
    if (block_a->starts[i] <= block_b->starts[i]) {
      intersect_block.starts[i] = block_b->starts[i];
      if (block_a->starts[i] + block_a->sizes[i] >
          block_b->starts[i] + block_b->sizes[i]) {
        intersect_block.sizes[i] = block_b->starts[i] + block_b->sizes[i] -
          intersect_block.starts[i];
      } else if (block_a->starts[i] + block_a->sizes[i] > block_b->starts[i]) {
        intersect_block.sizes[i] = block_a->starts[i] + block_a->sizes[i] -
          intersect_block.starts[i];
      } else {
        return 0;
      }
    } else {
      intersect_block.starts[i] = block_a->starts[i];
      if (block_b->starts[i] + block_b->sizes[i] >
          block_a->starts[i] + block_a->sizes[i]) {
        intersect_block.sizes[i] = block_a->starts[i] + block_a->sizes[i] -
          intersect_block.starts[i];
      } else if (block_b->starts[i] + block_b->sizes[i] > block_a->starts[i]) {
        intersect_block.sizes[i] = block_b->starts[i] + block_b->sizes[i] -
          intersect_block.starts[i];
      } else {
        return 0;
      }
    }
  }
  memcpy(ret_block, &intersect_block, sizeof(BIL_Block));
  ret_block->num_dims = block_a->num_dims;
  ret_block->total_size = 1;
  for (i = 0; i < ret_block->num_dims; i++) {
    ret_block->total_size *= ret_block->sizes[i];
  }
  return 1;
}

void BIL_Sched_compute_io_schedule(int num_dims,
                                   int* bb_mins, int* bb_maxs, int* file_res,
                                   int* var_sizes,
                                   int num_groups, char* group_names,
                                   int* blocks_to_group_map,
                                   BIL_Sched_IO* io_sched_ret,
                                   BIL_Sched_IO* inv_io_sched_ret) {
  int group_name_size = BIL_MAX_FILE_NAME_SIZE + BIL_MAX_VAR_NAME_SIZE;
  int i, j;
  io_sched_ret->num_io_blocks = 0;
  io_sched_ret->io_blocks = NULL;
  // Find how many io stages there will be.
  io_sched_ret->num_io_stages = ceil((float)num_groups / BIL->world_size);
  int s;
  int num_groups_to_assign = num_groups;
  int num_groups_assigned = 0;
  int* num_procs_per_group = 
    BIL_Misc_malloc(sizeof(int) * num_groups_to_assign);
  int* groups_root_proc =
    BIL_Misc_malloc(sizeof(int) * num_groups_to_assign);
  for (s = 0; s < io_sched_ret->num_io_stages; s++) {
    int num_procs_to_assign = BIL->world_size;
    num_groups = (num_groups_to_assign > num_procs_to_assign) ?
      num_procs_to_assign : num_groups_to_assign;
    num_groups_to_assign -= num_groups;

    // Find the amount of processes per group, the root process per group, and
    // your io group and rank.
    int io_group = -1;
    int io_group_rank = -1;
    int group_prefix_sum = 0;
    for (i = num_groups_assigned; i < num_groups_assigned + num_groups; i++) {
      num_procs_per_group[i] =
        num_procs_to_assign / (num_groups_assigned + num_groups - i);
      num_procs_to_assign -= num_procs_per_group[i];

      groups_root_proc[i] = group_prefix_sum; 
      group_prefix_sum += num_procs_per_group[i];
      if (BIL->world_rank < group_prefix_sum && io_group == -1) {
        io_group = i;
        io_group_rank = 
          BIL->world_rank - (group_prefix_sum - num_procs_per_group[i]);
      }
    }

    // Find the block of data that will be read for the I/O group
    BIL_Block io_block;
    BIL_Sched_get_io_block_assignment(
        num_dims, io_group, io_group_rank,
        num_procs_per_group[io_group],
        bb_mins + (io_group * BIL_MAX_NUM_DIMS),
        bb_maxs + (io_group * BIL_MAX_NUM_DIMS),
        file_res + (io_group * BIL_MAX_NUM_DIMS),
        group_names + (io_group * group_name_size),
        var_sizes[io_group], &io_block);
    // Compute misc information about the block.
    io_block.total_size = 1;
    for (j = 0; j < io_block.num_dims; j++) {
      io_block.total_size *= io_block.sizes[j];
    }

    // Add the io block to the io schedule.
    io_sched_ret->num_io_blocks++;
    io_sched_ret->io_blocks =
      BIL_Misc_realloc(io_sched_ret->io_blocks,
                       io_sched_ret->num_io_blocks * sizeof(BIL_Block));
    memcpy(io_sched_ret->io_blocks + io_sched_ret->num_io_blocks - 1,
           &io_block, sizeof(BIL_Block));

    num_groups_assigned += num_groups;
  }
  
  // Create an inverse I/O schedule. This schedule determines which processes
  // will read your blocks.
  inv_io_sched_ret->num_io_blocks = 0;
  inv_io_sched_ret->io_blocks = NULL;
  inv_io_sched_ret->num_io_stages = -1;
  for (i = 0; i < BIL->num_blocks; i++) {
    // Find the I/O group that your block belongs to.
    int io_group = blocks_to_group_map[i];
    int io_group_root_proc = groups_root_proc[io_group];

    BIL_Block* io_blocks =
      BIL_Misc_malloc(sizeof(BIL_Block) * num_procs_per_group[io_group]);
    BIL_Sched_get_io_block_assignment(
        num_dims, io_group, -1, 
        num_procs_per_group[io_group],
        bb_mins + (io_group * BIL_MAX_NUM_DIMS),
        bb_maxs + (io_group * BIL_MAX_NUM_DIMS),
        file_res + (io_group * BIL_MAX_NUM_DIMS),
        group_names + (io_group * group_name_size),
        var_sizes[io_group], io_blocks);
    
    for (j = 0; j < num_procs_per_group[io_group]; j++) {
      // Intersect your block with the computed I/O block
      BIL_Block inv_io_block;
      memcpy(&inv_io_block, &(BIL->blocks[i]), sizeof(BIL_Block));
      if (BIL_Sched_intersect_block(&(io_blocks[j]), &inv_io_block,
                                    &inv_io_block) != 0) {
        inv_io_block.read_rank = io_group_root_proc + j;
        inv_io_block.request_rank = BIL->world_rank;
        memcpy(inv_io_block.file_name, io_blocks[j].file_name,
               BIL_MAX_FILE_NAME_SIZE);
        memcpy(inv_io_block.var_name, io_blocks[j].var_name,
               BIL_MAX_VAR_NAME_SIZE);
        memcpy(inv_io_block.file_dim_sizes, io_blocks[j].sizes,
               sizeof(int) * BIL_MAX_NUM_DIMS);
        inv_io_block.block_id = i;
        inv_io_block.var_size = io_blocks[j].var_size;
        inv_io_block.data = NULL;

        // Add the io block to the io schedule.
        inv_io_sched_ret->num_io_blocks++;
        inv_io_sched_ret->io_blocks =
          BIL_Misc_realloc(inv_io_sched_ret->io_blocks,
                           inv_io_sched_ret->num_io_blocks
                             * sizeof(BIL_Block));
        memcpy(inv_io_sched_ret->io_blocks 
                 + inv_io_sched_ret->num_io_blocks - 1,
               &inv_io_block, sizeof(BIL_Block));
      }
    }
    
    BIL_Misc_free(io_blocks);
/*
    for (j = 0; j < num_procs_per_group[io_group]; j++) {
      BIL_Block io_block;
      BIL_Sched_get_io_block_assignment(
          num_dims, io_group, j, 
          num_procs_per_group[io_group],
          bb_mins + (io_group * BIL_MAX_NUM_DIMS),
          bb_maxs + (io_group * BIL_MAX_NUM_DIMS),
          file_res + (io_group * BIL_MAX_NUM_DIMS),
          group_names + (io_group * group_name_size),
          var_sizes[io_group], &io_block);
       
      // Intersect your block with the computed I/O block
      BIL_Block inv_io_block;
      memcpy(&inv_io_block, &(BIL->blocks[i]), sizeof(BIL_Block));
      if (BIL_Sched_intersect_block(&io_block, &inv_io_block, &inv_io_block)
          != 0) {
        inv_io_block.read_rank = io_group_root_proc + j;
        inv_io_block.request_rank = BIL->world_rank;
        memcpy(inv_io_block.file_name, io_block.file_name,
               BIL_MAX_FILE_NAME_SIZE);
        memcpy(inv_io_block.var_name, io_block.var_name,
               BIL_MAX_VAR_NAME_SIZE);
        memcpy(inv_io_block.file_dim_sizes, io_block.sizes,
               sizeof(int) * BIL_MAX_NUM_DIMS);
        inv_io_block.block_id = i;
        inv_io_block.var_size = io_block.var_size;
           
        // Add the io block to the io schedule.
        inv_io_sched_ret->num_io_blocks++;
        inv_io_sched_ret->io_blocks =
          BIL_Misc_realloc(inv_io_sched_ret->io_blocks,
                           inv_io_sched_ret->num_io_blocks
                             * sizeof(BIL_Block));
        memcpy(inv_io_sched_ret->io_blocks 
                 + inv_io_sched_ret->num_io_blocks - 1,
               &inv_io_block, sizeof(BIL_Block));
      }
    }
*/
  }
  BIL_Misc_free(num_procs_per_group);
  BIL_Misc_free(groups_root_proc);
}

void BIL_Sched_get(BIL_Sched_IO* io_sched_ret, BIL_Sched_IO* inv_io_sched_ret,
                   int* num_groups_ret, int** blocks_to_group_map_ret) {
  BIL_Timing_sched_start();
  
  // Gather all the unique group names. A group name consists of a file name
  // and a variable name.
  char* group_names;
  int num_groups;
  BIL_Sched_gather_group_names(&group_names, &num_groups);

  // After finding the unique group names, create a map from BIL's blocks to
  // these group names.
  int i;
  int* blocks_to_group_map = 
    BIL_Misc_malloc(sizeof(int) * BIL->num_blocks);
  int group_name_size = BIL_MAX_FILE_NAME_SIZE + BIL_MAX_VAR_NAME_SIZE;
  char* group_name = BIL_Misc_malloc(group_name_size);
  for (i = 0; i < BIL->num_blocks; i++) {
    memcpy(group_name, BIL->blocks[i].file_name, BIL_MAX_FILE_NAME_SIZE);
    memcpy(group_name + BIL_MAX_FILE_NAME_SIZE, BIL->blocks[i].var_name,
           BIL_MAX_VAR_NAME_SIZE);
    blocks_to_group_map[i] =
      ((uint64_t)bsearch(group_name, group_names, num_groups,
                         group_name_size,  BIL_Sched_compare_group_name) -
       (uint64_t)group_names) / group_name_size;
  }
  BIL_Misc_free(group_name);

  // Find out the bounding box and file resolution for each group.
  int* bb_mins = BIL_Misc_malloc(sizeof(int) * num_groups * BIL_MAX_NUM_DIMS);
  int* bb_maxs = BIL_Misc_malloc(sizeof(int) * num_groups * BIL_MAX_NUM_DIMS);
  int* file_res = BIL_Misc_malloc(sizeof(int) * num_groups * BIL_MAX_NUM_DIMS);
  memset(file_res, 0, sizeof(int) * num_groups * BIL_MAX_NUM_DIMS);
  int* var_size = BIL_Misc_malloc(sizeof(int) * num_groups);
  memset(var_size, 0, sizeof(int) * num_groups);
  for (i = 0; i < num_groups * BIL_MAX_NUM_DIMS; i++) {
    bb_mins[i] = INT_MAX;
    bb_maxs[i] = INT_MIN;
  }
  int k;
  int num_dims = -1;
  for (i = 0; i < BIL->num_blocks; i++) {
    int group = blocks_to_group_map[i];
    BIL->io_type = BIL->blocks[i].io_type;
    for (k = 0; k < BIL->blocks[i].num_dims; k++) {
      bb_mins[group * BIL_MAX_NUM_DIMS + k] =
        (BIL->blocks[i].starts[k] < bb_mins[group * BIL_MAX_NUM_DIMS + k]) ?
          BIL->blocks[i].starts[k] : bb_mins[group * BIL_MAX_NUM_DIMS + k];
      bb_maxs[group * BIL_MAX_NUM_DIMS + k] =
        (BIL->blocks[i].starts[k] + BIL->blocks[i].sizes[k] > 
         bb_maxs[group * BIL_MAX_NUM_DIMS + k]) ?
           BIL->blocks[i].starts[k] + BIL->blocks[i].sizes[k] :
           bb_maxs[group * BIL_MAX_NUM_DIMS + k];
      file_res[group * BIL_MAX_NUM_DIMS + k] = BIL->blocks[i].file_dim_sizes[k];
      num_dims = BIL->blocks[i].num_dims;
      var_size[group] = BIL->blocks[i].var_size;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &num_dims, 1, MPI_INT, MPI_MAX, BIL->world_comm);
  MPI_Allreduce(MPI_IN_PLACE, &(BIL->io_type), 1, MPI_INT, MPI_MAX,
                BIL->world_comm);
  MPI_Allreduce(MPI_IN_PLACE, bb_mins, num_groups * BIL_MAX_NUM_DIMS, MPI_INT,
                MPI_MIN, BIL->world_comm);
  MPI_Allreduce(MPI_IN_PLACE, bb_maxs, num_groups * BIL_MAX_NUM_DIMS, MPI_INT,
                MPI_MAX, BIL->world_comm);
  if (BIL->io_type == BIL_RAW) {
    MPI_Allreduce(MPI_IN_PLACE, file_res, num_groups * BIL_MAX_NUM_DIMS,
                  MPI_INT, MPI_MAX, BIL->world_comm);
    MPI_Allreduce(MPI_IN_PLACE, var_size, num_groups, MPI_INT, MPI_MAX,
                  BIL->world_comm);
  }

  // Compute the I/O schedule.
  BIL_Sched_compute_io_schedule(num_dims, bb_mins, bb_maxs, file_res, 
                                var_size, num_groups, group_names,
                                blocks_to_group_map, io_sched_ret,
                                inv_io_sched_ret);

  // Clean up.
  BIL_Misc_free(bb_mins);
  BIL_Misc_free(bb_maxs);
  BIL_Misc_free(file_res);
  BIL_Misc_free(var_size);
  BIL_Misc_free(group_names);

  *blocks_to_group_map_ret = blocks_to_group_map;
  *num_groups_ret = num_groups;

  BIL_Timing_sched_stop();
}
