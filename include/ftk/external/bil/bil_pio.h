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

// bil_pio.h
// Wes Kendall
// 07/18/2010
// Function prototypes for aggregating and issuing I/O in BIL.
////

#ifndef __BIL_PIO_H
#define __BIL_PIO_H

#include "bil.h"
#include "bil_sched.h"

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

void BIL_Pio_issue(BIL_Sched_IO* io_sched, BIL_Sched_IO* inv_io_sched,
                   int num_groups, int* blocks_to_group_map);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // __BIL_PIO_H
