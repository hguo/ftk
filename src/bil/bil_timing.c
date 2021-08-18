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

// bil_timing.c
// Wes Kendall
// 07/18/2010
// Functions for timing I/O, communication, and computation in BIL.
////

#include <stdio.h>

#include "bil.h"
#include "bil_timing.h"
#include "bil_misc.h"

void BIL_Timing_total_start() {
  MPI_Barrier(BIL->world_comm);
  BIL->timing.total_time -= MPI_Wtime();
}

void BIL_Timing_total_stop() {
  MPI_Barrier(BIL->world_comm);
  BIL->timing.total_time += MPI_Wtime();
}

void BIL_Timing_comm_start() {
  BIL->timing.comm_time -= MPI_Wtime();
}

void BIL_Timing_comm_stop() {
  BIL->timing.comm_time += MPI_Wtime();
}

void BIL_Timing_sched_start() {
  BIL->timing.sched_time -= MPI_Wtime();
}

void BIL_Timing_sched_stop() {
  BIL->timing.sched_time += MPI_Wtime();
}

void BIL_Timing_comp_start() {
  BIL->timing.comp_time -= MPI_Wtime();
}

void BIL_Timing_comp_stop() {
  BIL->timing.comp_time += MPI_Wtime();
}

void BIL_Timing_io_start(MPI_Comm timing_comm) {
#ifndef DISABLE_TIMING
  MPI_Barrier(timing_comm);
#endif
  BIL->timing.io_time -= MPI_Wtime();
}

void BIL_Timing_io_stop(MPI_Comm timing_comm, const int64_t bytes_read) {
#ifndef DISABLE_TIMING
  MPI_Barrier(timing_comm);
#endif
  BIL->timing.io_time += MPI_Wtime();
  BIL->timing.bytes_read += bytes_read;
}

void BIL_Timing_fopen_start(MPI_Comm timing_comm) {
#ifndef DISABLE_TIMING
  MPI_Barrier(timing_comm);
#endif
  BIL->timing.fopen_time -= MPI_Wtime();
}

void BIL_Timing_fopen_stop(MPI_Comm timing_comm) {
#ifndef DISABLE_TIMING
  MPI_Barrier(timing_comm);
#endif
  BIL->timing.fopen_time += MPI_Wtime();
}

// BIL_Timing_init
// Initializes timing values of BIL.
////
void BIL_Timing_init(BIL_Timing *timing) {
  timing->bytes_read = 0;
  timing->total_time = 0;
  timing->comp_time = 0;
  timing->comm_time = 0;
  timing->fopen_time = 0;
  timing->io_time = 0;
}

// BIL_Timing_print
// Prints timing statistics of BIL.
////
void BIL_Timing_print() {
  int64_t total_bytes_read;
  double max_io_time;
  MPI_Reduce(&(BIL->timing.io_time), &max_io_time, 1, MPI_DOUBLE, MPI_MAX,
             0, BIL->world_comm);
  double avg_io_time;
  MPI_Reduce(&(BIL->timing.io_time), &avg_io_time, 1, MPI_DOUBLE, MPI_SUM,
             0, BIL->world_comm);
  avg_io_time /= BIL->world_size;
  double max_open_time;
  MPI_Reduce(&(BIL->timing.fopen_time), &max_open_time, 1, MPI_DOUBLE, MPI_MAX,
             0, BIL->world_comm);
  double avg_open_time;
  MPI_Reduce(&(BIL->timing.fopen_time), &avg_open_time, 1, MPI_DOUBLE, MPI_SUM,
             0, BIL->world_comm);
  avg_open_time /= BIL->world_size;
  MPI_Reduce(&(BIL->timing.bytes_read), &total_bytes_read, 1, MPI_LONG_LONG,
      MPI_SUM, 0, BIL->world_comm);
  double max_comp_time;
  MPI_Reduce(&(BIL->timing.comp_time), &max_comp_time, 1, MPI_DOUBLE, MPI_MAX,
             0, BIL->world_comm);
  double avg_comp_time;
  MPI_Reduce(&(BIL->timing.comp_time), &avg_comp_time, 1, MPI_DOUBLE, MPI_SUM,
             0, BIL->world_comm);
  avg_comp_time /= BIL->world_size;
  double max_comm_time;
  MPI_Reduce(&(BIL->timing.comm_time), &max_comm_time, 1, MPI_DOUBLE, MPI_MAX,
             0, BIL->world_comm);
  double avg_comm_time;
  MPI_Reduce(&(BIL->timing.comm_time), &avg_comm_time, 1, MPI_DOUBLE, MPI_SUM,
             0, BIL->world_comm);
  avg_comm_time /= BIL->world_size;
  if (BIL->world_rank == 0) {
    if (max_io_time != 0.0) {
      fprintf(stderr, 
              "BIL_STATS:%d:MaxOpen:%lf:AvgOpen:%lf:MaxRead:%lf:AvgRead:%lf:"
              "Bytes:%lld:Bandwidth:%lf:MaxComm:%lf:AvgComm:%lf:MaxComp:%lf:"
              "AvgComp:%lf:Sched:%lf:Total:%lf:AggBandwidth:%lf\n",
              BIL->world_size, max_open_time, avg_open_time, max_io_time,
              avg_io_time, total_bytes_read,
              (total_bytes_read / (1024.0f * 1024.0f)) / max_io_time,
              max_comm_time, avg_comm_time, max_comp_time, avg_comp_time,
              BIL->timing.sched_time, BIL->timing.total_time,
              (total_bytes_read / (1024.0f * 1024.0f)) 
                / BIL->timing.total_time);
    } else {
      fprintf(stderr, "No timing information available from BIL\n");
    }
  }
}
