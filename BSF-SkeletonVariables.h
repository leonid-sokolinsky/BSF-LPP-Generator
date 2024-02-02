/*==============================================================================
Project: LiFe
Theme: LPP Generator
Module: BSF-SceletonVariables.h (Skeleton Variables)
Prefix: BSF_sv
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
static int BSF_sv_addressOffset;			// Contains the number of the first element of the map-sublist scheduled to the current worker process.
static int BSF_sv_iterCounter;				// Contains the number of iterations performed so far.
static int BSF_sv_jobCase;					// Contains the number of the current activity (job) in workflow.
static int BSF_sv_mpiRank;					// Contains the rank of current MPI process.
static int BSF_sv_mpiMaster;				// Contains the rank of master MPI process.
static int BSF_sv_numberInSublist;			// Contains the relative number of the element in the map-sublist that the function Map is currently applied to. 
static int BSF_sv_numOfWorkers;				// Contains the total number of the worker processes.
static PT_bsf_parameter_T BSF_sv_parameter; // Contains the order parameters.
static int BSF_sv_sublistLength;			// Contains the length of a map-sublist scheduled to a worker process.