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

// bil.h
// Wes Kendall
// 07/18/2010
// The prototypes for the main BIL interface and the only functions available
// to the user. A manual with descriptions about most of these functions is
// available at http://seelab.eecs.utk.edu/bil.php.
////

#ifndef __BIL_H
#define __BIL_H

#include <ftk/config.hh>
#if (!FTK_HAVE_PNETCDF)
  #ifndef DISABLE_PNETCDF
    #define DISABLE_PNETCDF
  #endif 
#endif

#include <mpi.h>
#include <inttypes.h>

#if HAVE_CONFIG_H
#include <config.h>
static const char version_number[] = VERSION;
#endif // HAVE_CONFIG_H

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

#define BIL_MAX_NUM_DIMS 4
#define BIL_MAX_FILE_NAME_SIZE 128
#define BIL_MAX_VAR_NAME_SIZE 128
#define BIL_TIMING

// BIL_Error_Type
// An enum for all of the internal errors of BIL.
////
typedef enum {
  BIL_BLOCK_ASSIGN_ERR = 0,
  BIL_ONE_PHASE_IO_ERR,
  BIL_BLOCK_DIM_MATCH_ERR,
  BIL_UNSUPPORTED_NC_TYPE_ERR,
  BIL_UNSUPPORTED_IO_TYPE_ERR,
  BIL_BLOCK_TOO_MANY_DIMS_ERR
} BIL_Error_Type;

// BIL_IO_Type
// An enum describing the types of I/O.
////
typedef enum {
  BIL_RAW = 0,
#ifndef DISABLE_PNETCDF
  BIL_NC,
#endif // DISABLE_PNETCDF
  BIL_NONE
} BIL_IO_Type;

// BIL_Timing 
// The main struct for timing information in BIL.
////
typedef struct {
  int64_t bytes_read; // Total bytes read.
  double total_time; // Total time.
  double io_time; // Time to perform I/O.
  double comm_time; // Time to perform communication.
  double comp_time; // Time to perform computation.
  double fopen_time; // Time to open files.
  double sched_time; // Time to perform scheduling.
} BIL_Timing;

// BIL_Block
// A structure for holding information about a block of data.
////
typedef struct {
  void* data; // Data buffer for reading/writing.
  void** user_data_buffer; // Pointer to the user's data buffer that they
                           // supplied. This is merely for convenience and to
                           // simplify the code since the user may either
                           // supply the buffer or have it allocated for them. 
  int num_dims; // Number of block dimensions.
  int starts[BIL_MAX_NUM_DIMS]; // Block dimension starts.
  int sizes[BIL_MAX_NUM_DIMS]; // Block dimension sizes.
  int total_size; // Total size of dimensions.
  int file_dim_sizes[BIL_MAX_NUM_DIMS]; // File dimension sizes.
  char file_name[BIL_MAX_FILE_NAME_SIZE]; // Filename.
  char var_name[BIL_MAX_VAR_NAME_SIZE]; // Variable name from file.
  int request_rank; // Rank of process that requested the I/O.
  int read_rank; // Rank of process that performs the read.
  int block_id; // Block id, local to the requesting process.
  int var_size; // Variable size.
  int io_group; // I/O group of process.
  BIL_IO_Type io_type; // Type of I/O for block.
} BIL_Block;

// BIL_Info
// Contains the global information used in BIL. This structure is initialized
// with BIL_Init and freed with BIL_Finalize.
////
typedef struct {
  MPI_Comm world_comm; // The world BIL is operating in.
  int world_rank; // World rank.
  int world_size; // World size.
  MPI_Datatype bil_block_type; // Datatype for BIL_Block struct.
  int num_blocks; // Num blocks requested by user.
  BIL_Block *blocks; // Blocks requested by user.
  BIL_Timing timing; // Additional timing information.
  BIL_IO_Type io_type; // The type of I/O that is being performed.
  // Parameters that the user may set.
  MPI_Info io_hints; // Hints to pass to MPI-I/O.
  int io_header_size; // Size of a raw header.
} BIL_Info;

// The main global structure of BIL.
extern BIL_Info *BIL;

// BIL_Init
// Initializes the BIL interface.
////
void BIL_Init(const MPI_Comm world_comm);

// BIL_Finalize
// Finalizes the BIL interface and frees all the structures.
////
void BIL_Finalize();

// BIL_Add_block_raw
// Adds a raw data block to the main global structure. A read is issued on this
// block when BIL_Read is called.
////
void BIL_Add_block_raw(int num_dims, const int *data_dims,
                       const int *block_start, const int *block_size,
                       const char *file_name, MPI_Datatype var_type,
                       void** buffer);

// BIL_Add_block_nc
// Adds a netcdf data block to the main global structure. A read is then issued
// on this block when BIL_Read is called.
////
#ifndef DISABLE_PNETCDF
void BIL_Add_block_nc(int num_dims, const int *block_start,
                      const int *block_size, const char *file_name,
                      const char *var_name, void** buffer);
#endif // DISABLE_PNETCDF

// BIL_Read
// The main function for reading data. This function is a wrapper around an I/O
// schedule generator and the actual I/O function.
////
void BIL_Read();

// BIL_Set_io_hints
// Sets the I/O hints that are used during parallel I/O.
////
void BIL_Set_io_hints(MPI_Info io_hints);

// BIL_Set_io_header_size();
// Sets the I/O header size of each file.
////
void BIL_Set_io_header_size(int io_header_size);

// BIL_Get_timing
// Gets the timing information for BIL.
////
void BIL_Get_timing();

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // __BIL_H
