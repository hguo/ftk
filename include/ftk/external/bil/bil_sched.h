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

// bil_sched.h
// Wes Kendall
// 07/18/2010
// Function prototypes for scheduling blocks to processes in BIL.
////

#ifndef __BIL_SCHED_H
#define __BIL_SCHED_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct {
  BIL_Block* io_blocks;  // The blocks being read.
  int num_io_blocks;  // The number of blocks to read.
  int is_two_phase;  // Is the schedule a two phase I/O schedule?
  BIL_Block* my_blocks_being_read; // My I/O blocks that are being read.
  int my_num_blocks_being_read;  // The number of my blocks that are being read.
} BIL_Sched;

typedef struct {
  BIL_Block* io_blocks;
  int num_io_blocks;
  int num_io_stages;
} BIL_Sched_IO;

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

int BIL_Sched_compare_file_and_var_name(const void* a, const void* b);

void BIL_Sched_finalize(BIL_Sched* bil_sched);

void BIL_Sched_get(BIL_Sched_IO* io_sched_ret, BIL_Sched_IO* inv_io_sched_ret,
                   int* num_groups_ret, int** blocks_to_group_map_ret);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // __BIL_SCHED_H
