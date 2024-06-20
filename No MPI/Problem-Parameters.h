/*==============================================================================
Project: LiFe
Theme: LPP Generator (no MPI)
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author: Leonid B. Sokolinsky
Publication: Sokolinsky L.B., Sokolinskaya I.M. FRaGenLP: A Generator of Random Linear 
	Programming Problems for Cluster Computing Systems // Parallel Computational Technologies. 
	PCT 2021. Communications in Computer and Information Science. 2021, vol. 1437. 164-177. 
	DOI:10.1007/978-3-030-81691-9_12.
==============================================================================*/
//#define PP_DEBUG
#define PP_PATH "D:/YandexDisk/_private/Programming/Set-of-LP-Problems/Rnd-LP/"
//=========================== Problem Parameters =========================
#define PP_N 15
#define PP_NUM_OF_RND_INEQUALITIES 1	// Number of random inequalities		
#define PP_RND_SEED 6					// Value used by srand() to seed pseudo-random number generator rand(). 
										//		0 corresponds to value depending on time.
#define PP_M (PP_N + 1 + PP_NUM_OF_RND_INEQUALITIES) // Total number of inequalities (+-1)
#define PP_MTX_NON (PP_M * (PP_N + 1))	// Number of non-zero elements in matrix A for MTX format
#define PP_ALPHA 200					// Length of hypercube edge
#define PP_THETA (PP_ALPHA/2)			// Radius of large hypersphere
#define PP_RHO (PP_THETA/2)				// Radius of small hypersphere
#define	PP_A_MAX (PP_RHO/2)				// Maximal random value for A
#define	PP_LIKE_FACTOR 0.4				// Range of values [0, 0.7]. Lower value implies greater likeness of hyperplanes
#define	PP_MIN_SHIFT (PP_RHO/3)			// Minimal acceptable shift between hyperplanes
//-------------------------- Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	8	// Number of Elements to output
#define PP_SETW 8
//#define PP_MATRIX_OUTPUT	
//#define PP_PATH "D:/YandexDisk/_private/Programming/LP-Problems/"
//#define PP_PATH ""
#define PP_PROBLEM_NAME "rnd"
//------------------------- Matrix format ----------------
#define PP_MTX_PREFIX			"lp_"
#define PP_MTX_POSTFIX_A		".mtx"
#define PP_MTX_POSTFIX_B		"_b.mtx"
#define PP_MTX_POSTFIX_LO		"_lo.mtx"
#define PP_MTX_POSTFIX_HI		"_hi.mtx"
#define PP_MTX_POSTFIX_C		"_c.mtx"