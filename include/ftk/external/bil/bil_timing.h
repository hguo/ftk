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

// bil_timing.h
// Wes Kendall
// 07/18/2010
// Function prototypes for timing I/O, communication, and computation in BIL.
////

#ifndef __BIL_TIMING_H
#define __BIL_TIMING_H

#include "bil.h"
#include <inttypes.h>

#include <inttypes.h>
#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

void BIL_Timing_init(BIL_Timing *timing);

void BIL_Timing_print();

void BIL_Timing_total_start();

void BIL_Timing_total_stop();

void BIL_Timing_comm_start();

void BIL_Timing_comm_stop();

void BIL_Timing_sched_start();

void BIL_Timing_sched_stop();

void BIL_Timing_comp_start();

void BIL_Timing_comp_stop();

void BIL_Timing_io_start(MPI_Comm timing_comm);

void BIL_Timing_io_stop(MPI_Comm timing_comm, const int64_t bytes_read);

void BIL_Timing_fopen_start(MPI_Comm timing_comm);

void BIL_Timing_fopen_stop(MPI_Comm timing_comm);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // __BIL_TIMING_H
