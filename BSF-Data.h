/*==============================================================================
Project: LiFe
Theme: LPP Generator
Module: BSF-Data.h (Problem Independent Variables & Data Structures)
Prefix: BD
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "BSF-Include.h"			// Problem Independent "Include" Files
#include "BSF-Types.h"				// Problem Independent Types
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
// ========================= Constants ========================
#define BD_EXIT	true			// End of Work

//======================= Compatibility ==========================
#ifndef PP_BSF_MAX_MPI_SIZE
#ifdef PP_MAX_MPI_SIZE
#define PP_BSF_MAX_MPI_SIZE PP_MAX_MPI_SIZE
#endif // PP_MAX_MPI_SIZE
#endif // PP_BSF_MAX_MPI_SIZE

// ========================= Variables ========================
static int	BD_elemsPerWorker;	// Number of elements per worker
static bool BD_exit;			// Indicator of ending the calculations
static int	BD_jobCase = 0;		// Default value
static int	BD_listSize;		// Size of the list
static int	BD_masterRank;		// Rank of Master = Number of MPI processes - 1
static int	BD_newJobCase = 0;	// Default value
static int	BD_numOfWorkers;	// Number of Workers = Number of MPI processes - 1
static int	BD_rank;			// Rank of the current process MPI 
static int	BD_size;			// Number of MPI processes
static bool BD_success;			// Indicator of Successfulness of Initialization
static int	BD_tailLength;		// Length of the remainder of the list after dividing by the number of workers
// ========================= Time Variables ========================
static double BD_t;				// Total time
static int BD_iterCounter;
// ========================= Lists ========================
static BT_extendedReduceElem_T*		BD_extendedReduceList;
static BT_extendedReduceElem_T_1*	BD_extendedReduceList_1;
static BT_extendedReduceElem_T_2*	BD_extendedReduceList_2;
static BT_extendedReduceElem_T_3*	BD_extendedReduceList_3;
static BT_extendedReduceElem_T*		BD_extendedReduceResult_P;
static BT_extendedReduceElem_T_1*	BD_extendedReduceResult_P_1;
static BT_extendedReduceElem_T_2*	BD_extendedReduceResult_P_2;
static BT_extendedReduceElem_T_3*	BD_extendedReduceResult_P_3;
static PT_bsf_mapElem_T*			BD_mapSubList;
static BT_order_T					BD_order;		// Order for Workers
static MPI_Status					BD_status[PP_BSF_MAX_MPI_SIZE];		// MPI ststus
static MPI_Request					BD_request[PP_BSF_MAX_MPI_SIZE];	// MPI request
static int							BD_sublistSize[PP_BSF_MAX_MPI_SIZE];// SubList Sizes for Workers
static int							BD_offset[PP_BSF_MAX_MPI_SIZE];		// List offsets for Workers