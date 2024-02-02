/*==============================================================================
Project: LiFe
Theme: LPP Generator (no MPI)
Module: Problem-Data.h (Problem Data)
Prefix: PD
Author: Leonid B. Sokolinsky
==============================================================================*/
#include "Problem-Types.h"			// Problem Parameters 
using namespace std;

//========================== Problem variables ====================================
static int PD_n; // Current dimension
static int PD_m; // Current number of inequalities
static int PD_m_predef; // Number of predefined inequalities
static PT_float_T PD_sqrt_n; // Square root of n
static int PD_k; // Index of current random inequality
static PT_float_T PD_centerObjectF;	// Value of object function in the center of hypercube
static unsigned PD_failuresType1; // Hyperplane is similar to another one

//========================== Problem data structures ==============================
static PT_matrix_T PD_A;
static PT_column_T PD_b;
static PT_vector_T PD_c;
static PT_vector_T PD_center;		// Center of hypercube
static PT_column_T PD_aNorm;
static PT_MTX_A_nor_noc PD_MTX_A_nor_noc;
static PT_MTX_A_val PD_MTX_A_val;

//========================== Files ==============================
static string PD_MTX_File;
static string PD_problemName;