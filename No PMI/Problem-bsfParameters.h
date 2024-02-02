/*==============================================================================
Project: LiFe
Theme: LPP Generator (no MPI)
Module: Problem-bsfParameters.h (BSF-skeleton parameters)
Prefix: PP_BSF
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/

//=========================== Skeleton Parameters =========================
#define PP_BSF_MAX_MPI_SIZE 2		// Do not modify!
//#define PP_BSF_ITER_OUTPUT			// If it is defined then Iteration Output is performed
#define PP_BSF_TRACE_COUNT 1		// Each PP_BSF_TRACE_COUNT-th iteration to be outputted
#define PP_BSF_MAX_JOB_CASE 0		// Defines the maximum number of activities (jobs) in workflow minus 1
//--------------------------- OpenMP Parameters ---------------------------
#define PP_BSF_OMP				// If PP_BSF_OMP is defined then OpenMP is turned on for Map Step
//#define PP_BSF_NUM_THREADS 12		// If PP_BSF_NUM_THREADS is udefined then all accessable threads are used
//--------------- BSF Lists parameters (For "No MPI" only) ----------------
#include "Problem-Parameters.h"
#define PP_BSF_MAP_LIST_LENGTH		1
#define PP_BSF_REDUCE_LIST_LENGTH	1
#define PP_BSF_REDUCE_LIST_1_LENGTH	1
#define PP_BSF_REDUCE_LIST_2_LENGTH	1
#define PP_BSF_REDUCE_LIST_3_LENGTH	1