/*==============================================================================
Project: LiFe
Theme: LPP Generator (no MPI)
Module: Problem-bsfCode.cpp (Problem-dependent Code)
Prefix: PC
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;

//----------------------- Predefined problem-dependent functions -----------------
void PC_bsf_Init(bool* success) {

	if (PP_BSF_MAX_MPI_SIZE != 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "PP_BSF_MAX_MPI_SIZE must be equal to 2!" << endl;
		*success = false;
		return;
	}

	if (PP_RND_SEED > 0)
		srand(PP_RND_SEED);
	else
		srand((unsigned)time(NULL) * (BSF_sv_mpiRank + 10));

	PD_n = PP_N;

	if (PD_n < 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "PP_N must be greater than 1!" << endl;
		*success = false;
		return;
	}

	PD_m = PD_n;
	PD_k = PD_n;

	for (int i = 0; i < PD_m; i++) {
		for (int j = 0; j < PD_n; j++)
			PD_A[i][j] = 0;
		PD_A[i][i] = 1;
		PD_b[i] = PP_ALPHA;
	}

	if (PP_NUM_OF_RND_INEQUALITIES == 0) { // Standart problem
		for (int j = 0; j < PD_n; j++)
			PD_A[PD_m][j] = 1;
		PD_b[PD_m] = PP_ALPHA * (PD_m - 1) + (PT_float_T)PP_ALPHA / 2;
		PD_m++; assert(PD_m <= PP_M);
		PD_k++; assert(PD_k <= PP_M);
	}

	PD_m_predef = PD_m;

	for (int j = 0; j < PD_n; j++)
		PD_center[j] = PP_ALPHA / 2;

	if (PP_NUM_OF_RND_INEQUALITIES == 0) // Standart objective function
		for (int j = 0; j < PD_n; j++)
			PD_c[j] = (PT_float_T)(j + 1);
	else // Random objective function
	{
		for (int j = 0; j < PD_n; j++)
			PD_c[j] = 0;
		for (int j = 0; j < PD_n; j++) {
			int rnd_j = rand() % PD_n;
			while (PD_c[rnd_j] != 0)
				rnd_j = rand() % PD_n;
			PD_c[rnd_j] = (PT_float_T)(j + 1);
		}
	}

	for (int i = 0; i < PD_m; i++)
		PD_aNorm[i] = sqrt(Vector_NormSquare(PD_A[i]));

	PD_sqrt_n = (PT_float_T)sqrt(PD_n);

	PD_problemName = PP_PROBLEM_NAME;
	
	PD_problemName += to_string(PD_n);
	PD_problemName += "-";
	PD_problemName += to_string(PP_NUM_OF_RND_INEQUALITIES);
	if (PP_RND_SEED > 0) {
		PD_problemName += "-";
		PD_problemName += to_string(PP_RND_SEED);
	}
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = BSF_sv_numOfWorkers;
}

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {

}

void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success) {	// For Job 0
	PT_float_T aNormSquare;

	reduceElem->failuresType1 = 0; // Hyperplane is similar to another one
	reduceElem->failuresType2 = 0; // Point (0,...,0,PP_ALPHA,0,...,0) is not feasible

	if (PP_NUM_OF_RND_INEQUALITIES == 0)
		return;

	do {
		for (int j = 0; j < PD_n; j++) // computing a[*]
			reduceElem->a[j] = PP_THETA + RndSign() * RndValue(PP_THETA);

		PT_float_T term; // computing b
		term = PP_RHO / PD_sqrt_n + PP_ALPHA / 2 + RndValue((PP_THETA - PP_RHO) / PD_sqrt_n);
		reduceElem->b = 0;
		for (int j = 0; j < PD_n; j++)
			reduceElem->b += reduceElem->a[j] * term;

		bool failure = false;
		for (int j = 0; j < PD_n; j++)
			if (reduceElem->b / reduceElem->a[j] <= PP_ALPHA) { // Point (0,...,0,PP_ALPHA,0,...,0) is not feasible
				failure = true;
				break;
			}
		if (failure) {
			reduceElem->failuresType2++;
			continue;
		}

		aNormSquare = Vector_NormSquare(reduceElem->a);
		reduceElem->aNorm = sqrt(aNormSquare);

		bool like = false;
		for (int i = 0; i < PD_m_predef; i++)
			if (like = Like(reduceElem->a, reduceElem->b, reduceElem->aNorm, PD_A[i], PD_b[i], PD_aNorm[i]))
				break;

		if (like) {// Hyperplane is similar to another one
			reduceElem->failuresType1++;
			continue;
		}
		
		break;
	} while (true);
}

void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem, int* success) {// For Job 1
	// Not used
}

void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem, int* success) {// For Job 2
	// Not used
}

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {// For Job 3
	// Not used
}

void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) {			// For Job 0
	z->failuresType1 = x->failuresType1 + y->failuresType1;
	z->failuresType2 = x->failuresType2 + y->failuresType2;
}

void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {	// For Job 1
	// Not used
}

void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {	// For Job 2
	// Not used
}

void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {	// For Job 3
	// Not used
}

void PC_bsf_ProcessResults(		// For Job 0
	PT_bsf_reduceElem_T* reduceResult,
	int reduceCounter, 
	PT_bsf_parameter_T* parameter, 
	int* nextJob,
	bool* exit 
) {
	bool like;
	struct extendedReduceElem_T {	// Extended element type of reduce list
		PT_bsf_reduceElem_T elem;	// Element of reduce list
		int reduceCounter;			// Reduce Counter
	};
	extendedReduceElem_T* extendedReduceElem;

	if (PD_k == PP_NUM_OF_RND_INEQUALITIES + PD_m_predef) {
		*exit = true;
		return;
	}

	PD_failuresType1 += reduceResult->failuresType1;
	PD_failuresType2 += reduceResult->failuresType2;

	extendedReduceElem = (extendedReduceElem_T*)reduceResult;

	for (int w = 0; w < BSF_sv_numOfWorkers; w++) {

		like = false;

		for (int i = PD_n + 1; i < PD_k; i++) {
			if (like = Like(extendedReduceElem[w].elem.a, extendedReduceElem[w].elem.b, extendedReduceElem[w].elem.aNorm, PD_A[i], PD_b[i], PD_aNorm[i]))
				break;
		}

		if (like) {
			PD_failuresType1++;
			continue;
		}

		Vector_Copy(extendedReduceElem[w].elem.a, PD_A[PD_k]);
		PD_b[PD_k] = extendedReduceElem[w].elem.b;
		PD_aNorm[PD_k] = extendedReduceElem[w].elem.aNorm;
		PD_k++; assert(PD_k <= PP_M);

		if (PD_k == PP_NUM_OF_RND_INEQUALITIES + 2 * PD_n + 1) {
			*exit = true;
			return;
		}
	}
}

void PC_bsf_ProcessResults_1(	// For Job 1	
	PT_bsf_reduceElem_T_1* reduceResult,
	int reduceCounter, 
	PT_bsf_parameter_T* parameter, 
	int* nextJob,
	bool* exit 
) {
	// Not used
}

void PC_bsf_ProcessResults_2(	// For Job 2
	PT_bsf_reduceElem_T_2* reduceResult,
	int reduceCounter, 
	PT_bsf_parameter_T* parameter, 
	int* nextJob,
	bool* exit 
	) {
	// Not used
}

void PC_bsf_ProcessResults_3(	// For Job 3
	PT_bsf_reduceElem_T_3* reduceResult,
	int reduceCounter, 
	PT_bsf_parameter_T* parameter, 
	int* nextJob,
	bool* exit 
	) {
	// Not used
}

void PC_bsf_JobDispatcher(
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* job,
	bool* exit
) {
	// Not used
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "-------------------------------------PC_bsf_ParametersOutput-----------------------------------" << endl;
	cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
#else
	cout << "Number of Threads: " << omp_get_num_procs() << endl;
#endif // PP_BSF_NUM_THREADS
#else
	cout << "OpenMP is turned off!" << endl;
#endif // PP_BSF_OMP
	cout << "Dimension n = " << PP_N << endl;
	cout << "Number of random inequalities: " << PP_NUM_OF_RND_INEQUALITIES << endl;
	cout << "Length of hypercube edge ALPHA = " << PP_ALPHA << endl;
	cout << "Radius of large hypersphere RHO = " << PP_THETA << endl;
	cout << "Radius of small hypersphere THETA = " << PP_RHO << endl;
	cout << "Maximal acceptable likeness of equations MAX_LIKE = " << PP_LIKE_FACTOR << endl;
	cout << "Minimal acceptable shift MIN_SHIFT = " << PP_MIN_SHIFT << endl;
#ifdef PP_MATRIX_OUTPUT
	cout << "------- Support inequalities -------" << endl;
	for (int i = 0; i < PD_m; i++) {
		cout << i << ")";
		for (int j = 0; j < PD_n; j++)
			cout << setw(PP_SETW) << PD_A[i][j] << "\t";
		cout << "\t<=\t" << setw(PP_SETW) << PD_b[i] << endl;
	}
	cout << "n = " << PD_n << "\tm = " << PD_m << "\tk = " << PD_k << endl;
#endif // PP_MATRIX_OUTPUT

	cout << "Objective Function:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(PP_SETW) << PD_c[j];
	cout << (PP_OUTPUT_LIMIT < PP_N ? "	..." : "") << endl;
}

void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {	// For Job 0
	static int k = PD_m;
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;

	if (PD_k > k) {
		cout << PD_k - 1 << ")\t";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_A[PD_k - 1][j] << "\t";
		cout << (PP_OUTPUT_LIMIT < PD_n ? "	..." : "") << "<=\t" << setw(PP_SETW) << PD_b[PD_k - 1] << endl;
		k = PD_k;
	}
	cout << "Failures 'Similar' = " << PD_failuresType1 << endl;
	cout << "Failures '(0," << PP_ALPHA << ",0)' = " << PD_failuresType2 << endl;
	cout << "-------------------------------------" << endl;
}

void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {	// For Job 1
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	// Not used
}

void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {	// For Job 2
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	// Not used
}

void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {	// For Job 3
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	// Not used
}

void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) 
{
	PD_m = PD_k;
	cout << "=============================================" << endl;
	cout << "Time: " << t << endl;
	cout << "Iterations: " << BSF_sv_iterCounter << endl;

	if (PP_NUM_OF_RND_INEQUALITIES > 0) {
#ifdef PP_MATRIX_OUTPUT
		cout << "------- Random inequalities -------" << endl;
		for (int i = PD_m_predef; i < PD_m; i++) {
			cout << i << ")\t";
			for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_A[i][j] << "\t";
			cout << (PP_OUTPUT_LIMIT < PD_n ? "	..." : "") << "<=\t" << setw(PP_SETW) << PD_b[i] << endl;
		}
		cout << "-----------------------------------" << endl;
#endif // PP_MATRIX_OUTPUT
		cout << "Failures 'Similar' = " << PD_failuresType1 << endl;
		cout << "Failures '(0," << PP_ALPHA << ",0)' = " << PD_failuresType2 << endl;
	}

	MTX_AddFreeVariables();

#ifdef PP_DEBUG
#ifdef PP_MATRIX_OUTPUT
	cout << "-------------------------------------PC_bsf_ProblemOutput-----------------------------------" << endl;
	for (int i = 0; i < PD_m; i++) {
		cout << i << ")";
		for (int j = 0; j < PD_n; j++)
			cout << setw(PP_SETW) << PD_A[i][j] << "\t";
		cout << "\t<=\t" << setw(PP_SETW) << PD_b[i] << endl;
	}
	cout << "n = " << PD_n << "\tm = " << PD_m << "\tk = " << PD_k << endl;
#endif // PP_MATRIX_OUTPUT
#endif // PP_DEBUG

	cout << "-----------------------------------" << endl;

	const char* fileName;
	FILE* stream;

	// Creating lp_<PD_problemName>.mtx
	int non; // number of non-zero elements
	MTX_Make_A(PD_MTX_A_nor_noc, PD_MTX_A_val, &non);

	PD_MTX_File = PP_PATH;
	PD_MTX_File += PP_MTX_PREFIX;
	PD_MTX_File += PD_problemName;
	PD_MTX_File += PP_MTX_POSTFIX_A;
	fileName = PD_MTX_File.c_str();
	stream;
	stream = fopen(fileName, "w");
	if (stream == NULL) {
		cout << "Failure of opening file " << fileName << "!\n";
		return;
	}
	fprintf(stream, "%d %d %d\n", PD_m, PD_n, non);
	for (int k = 0; k < non; k++) 
		fprintf(stream, "%d %d %.16f\n", PD_MTX_A_nor_noc[k][0], PD_MTX_A_nor_noc[k][1], PD_MTX_A_val[k]);
	fclose(stream);
	cout << "File " << fileName << " successfully created." << endl;

	// Creating lp_<PD_problemName>_b.mtx
	PD_MTX_File = PP_PATH;
	PD_MTX_File += PP_MTX_PREFIX;
	PD_MTX_File += PD_problemName;
	PD_MTX_File += PP_MTX_POSTFIX_B;
	fileName = PD_MTX_File.c_str();
	stream = fopen(fileName, "w");
	if (stream == NULL) {
		cout << "Failure of opening file " << fileName << "!\n";
		return;
	}
	fprintf(stream, "%d 1\n", PD_m);
	for (int i = 0; i < PD_m; i++)
		fprintf(stream, "%.16f\n", PD_b[i]);
	fclose(stream);
	cout << "File " << fileName << " successfully created." << endl;

	// Creating lp_<PD_problemName>_c.mtx
	PD_MTX_File = PP_PATH;
	PD_MTX_File += PP_MTX_PREFIX;
	PD_MTX_File += PD_problemName;
	PD_MTX_File += PP_MTX_POSTFIX_C;
	fileName = PD_MTX_File.c_str();
	stream = fopen(fileName, "w");
	if (stream == NULL) {
		cout << "Failure of opening file " << fileName << "!\n";
		return;
	}
	fprintf(stream, "%d 1\n", PD_n);
	for (int j = 0; j < PP_N; j++)
		fprintf(stream, "%f\n", -PD_c[j]);
	for (int j = PP_N; j < PD_n; j++)
		fprintf(stream, "0\n");

	fclose(stream);
	cout << "File " << fileName << " successfully created." << endl;

	// Creating lp_<PD_problemName>_hi.mtx
	PD_MTX_File = PP_PATH;
	PD_MTX_File += PP_MTX_PREFIX;
	PD_MTX_File += PD_problemName;
	PD_MTX_File += PP_MTX_POSTFIX_HI;
	fileName = PD_MTX_File.c_str();
	stream = fopen(fileName, "w");
	if (stream == NULL) {
		cout << "Failure of opening file " << fileName << "!\n";
		return;
	}
	fprintf(stream, "%d 1\n", PD_n);
	for (int j = 0; j < PD_n; j++)
		fprintf(stream, "1e+308\n");
	fclose(stream);
	cout << "File " << fileName << " successfully created." << endl;

	// Creating lp_<PD_problemName>_lo.mtx
	PD_MTX_File = PP_PATH;
	PD_MTX_File += PP_MTX_PREFIX;
	PD_MTX_File += PD_problemName;
	PD_MTX_File += PP_MTX_POSTFIX_LO;
	fileName = PD_MTX_File.c_str();
	stream = fopen(fileName, "w");
	if (stream == NULL) {
		cout << "Failure of opening file " << fileName << "!\n";
		return;
	}
	fprintf(stream, "%d 1\n", PD_n);
	for (int j = 0; j < PD_n; j++)
		fprintf(stream, "0\n");
	fclose(stream);
	cout << "File " << fileName << " successfully created." << endl;

	cout << "-----------------------------------" << endl;
}

void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {	// For Job 1
	// Not used
}

void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {	// For Job 2
	// Not used
}

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {	// For Job 3
	// Not used
}

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {

}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {

}

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; }
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; }
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; }
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; }
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; }
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; }
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; }
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; }

//----------------------------- User functions -----------------------------
inline PT_float_T Vector_NormSquare(PT_vector_T x) {
	PT_float_T s = 0;

	for (int j = 0; j < PD_n; j++)
		s += x[j] * x[j];
	return s;
}

inline void RndVector(PT_vector_T vector) {
	for (int j = 0; j < PD_n; j++)
		vector[j] = RndValue(PP_A_MAX);
}

inline PT_float_T RndValue(PT_float_T rndMax) { // rnd >= 0
	return ((PT_float_T)rand() / ((PT_float_T)RAND_MAX + 1)) * rndMax;
}

inline int RndSign() {
	int res = rand() % 2;
	if (res == 0)
		res = -1;
	return res;
}

inline bool Like(PT_vector_T a1, PT_float_T b1, PT_float_T a1Norm, PT_vector_T a2, PT_float_T b2, PT_float_T a2Norm) {
	PT_float_T like, shift;
	PT_vector_T e1, e2, e1_e2;

	Vector_MultiplyByNumber(a1, 1 / a1Norm, e1);
	Vector_MultiplyByNumber(a2, 1 / a2Norm, e2);

	Vector_Subtraction(e1, e2, e1_e2);
	like = sqrt(Vector_NormSquare(e1_e2));
	shift = fabs(b1 / a1Norm - b2 / a2Norm);
	if (like < PP_LIKE_FACTOR)
		if (shift < PP_MIN_SHIFT)
			return true;
	return false;
}

inline void Vector_MultiplyByNumber(PT_vector_T x, PT_float_T r, PT_vector_T y) {  // y = r*x
	for (int j = 0; j < PD_n; j++) {
		y[j] = x[j] * r;
	}
}

inline void Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
	for (int j = 0; j < PD_n; j++) {
		z[j] = x[j] - y[j];
	}
}

inline void Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
	for (int j = 0; j < PD_n; j++) {
		toPoint[j] = fromPoint[j];
	}
}

inline void MTX_AddFreeVariables() {
	for (int i = 0; i < PD_m; i++) {
		PD_A[i][i + PD_n] = 1;
	}
	PD_n += PD_m;
}

inline void MTX_Make_A(PT_MTX_A_nor_noc A_nor_noc, PT_MTX_A_val A_val, int* non) {
		*non = 0;
	for (int i = 0; i < PD_m; i++)
		for (int j = 0; j < PD_n; j++)
			if (PD_A[i][j] != 0) {
				A_nor_noc[*non][0] = i + 1;
				A_nor_noc[*non][1] = j + 1;
				A_val[*non] = PD_A[i][j];
				(*non)++; assert(*non < PP_MTX_NON);
			}
}