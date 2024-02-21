/*==============================================================================
Project: LiFe
Theme: LPP Generator (MPI)
Module: Problem-bsf-Forwards.h (Problem Function Forwards)
Author: Leonid B. Sokolinsky 
==============================================================================*/
#include "Problem-bsfTypes.h"	// Predefined BSF types
#include "Problem-Types.h"		// Problem Types
//====================== Problem Functions ===========================
bool		Like(PT_vector_T a1, PT_float_T b1, PT_float_T a1Norm, PT_vector_T a2, PT_float_T b2, PT_float_T a2Norm);
void		MTX_AddFreeVariables();
void		MTX_Make_A(PT_MTX_A_nor_noc A_nor_noc, PT_MTX_A_val A_val, int* non);
int			RndSign();
PT_float_T	RndValue(PT_float_T rndMax);
void		RndVector(PT_vector_T vector);
void		Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint);
void		Vector_MultiplyByNumber(PT_vector_T x, PT_float_T r, PT_vector_T y);
PT_float_T	Vector_NormSquare(PT_vector_T x);
void		Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z);
//====================== Macros ================================
#define PF_MIN(x,y) (x<y?x:y)
#define PF_MAX(x,y) (x>y?x:y)