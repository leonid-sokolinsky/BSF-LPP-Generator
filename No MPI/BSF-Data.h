/*==============================================================================
Project: LiFe
Theme: LPP Generator (no MPI)
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

/*No MPI*/
/*No MPI*/
/*No MPI*/
/*No MPI*/
/*No MPI*/
/*No MPI*/

// ========================= Variables ========================
static int	BD_elemsPerWorker;	// Number of elements per worker
static bool BD_exit;			// Indicator of ending the calculations
static int	BD_jobCase = 0;		// Default value
static int	BD_listSize;		// Size of the list
static int	BD_masterRank;		// Rank of Master = Number of MPI processes - 1
static int	BD_newJobCase = 0;	// Default value
static int	BD_numOfWorkers;	// Number of Workers = Number of MPI processes - 1
/*No MPI*/
/*No MPI*/
static bool BD_success;			// Indicator of Successfulness of Initialization
static int	BD_tailLength;		// Length of the remainder of the list after dividing by the number of workers
// ========================= Time Variables ========================
static double BD_t;				// Total time
static int BD_iterCounter;
// ========================= Lists ========================
static BT_extendedReduceElem_T		BD_extendedReduceList[PP_BSF_REDUCE_LIST_LENGTH];
static BT_extendedReduceElem_T_1	BD_extendedReduceList_1[PP_BSF_REDUCE_LIST_1_LENGTH];
static BT_extendedReduceElem_T_2	BD_extendedReduceList_2[PP_BSF_REDUCE_LIST_2_LENGTH];
static BT_extendedReduceElem_T_3	BD_extendedReduceList_3[PP_BSF_REDUCE_LIST_3_LENGTH];
static BT_extendedReduceElem_T*		BD_extendedReduceResult_P;
static BT_extendedReduceElem_T_1*	BD_extendedReduceResult_P_1;
static BT_extendedReduceElem_T_2*	BD_extendedReduceResult_P_2;
static BT_extendedReduceElem_T_3*	BD_extendedReduceResult_P_3;
static PT_bsf_mapElem_T				BD_mapList[PP_BSF_MAP_LIST_LENGTH];
static BT_order_T					BD_order;		// Order for Workers