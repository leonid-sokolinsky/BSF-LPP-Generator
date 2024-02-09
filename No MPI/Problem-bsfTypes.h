/*==============================================================================
Project: LiFe
Theme: LPP Generator (no MPI)
Module: Problem-bsfTypes.h (Predefined Problem-depended BSF Types)
Prefix: PT_bsf
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#pragma once
#include "Problem-Types.h"		// Problem Types 
//=========================== BSF Types =========================
struct PT_bsf_parameter_T {		// Order parameters
};

struct PT_bsf_mapElem_T {		// Type of map-list elements
};

struct PT_bsf_reduceElem_T {	// Type of reduce-list elements for Job 0 (default)	
	PT_vector_T a;
	PT_float_T b;
	PT_float_T aNorm;
	unsigned failuresType1; // Hyperplane is similar to another one
	unsigned failuresType2; // Point (0,...,0,PP_ALPHA,0,...,0) is not feasible
};

struct PT_bsf_reduceElem_T_1 {	// Type of reduce-list elements for Job 1
	// Not used
};

struct PT_bsf_reduceElem_T_2 {	// Type of reduce-list elements for Job 2
	// Not used
};

struct PT_bsf_reduceElem_T_3 {	// Type of reduce-list elements for Job 3
	// Not used
};