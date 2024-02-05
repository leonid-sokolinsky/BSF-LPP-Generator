/*==============================================================================
Project: LiFe
Theme: LPP Generator (no MPI)
Module: BSF-Code.cpp (Problem Independent Code)
Prefix: BC
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "BSF-Data.h"				// Problem Independent Variables & Data Structures 
#include "BSF-Forwards.h"			// Problem Independent Function Forwards
#include "BSF-ProblemFunctions.h"	// Predefined Problem Function Forwards
using namespace std;
//======================================== Problem-independent codes (don't modify them) ====================================
int main(int argc, char* argv[]) {
	PC_bsfAssignMpiRank(0);
	PC_bsfAssignMpiMaster(0);
	PC_bsfAssignNumOfWorkers(1);

	BD_success = true;
	PC_bsf_Init(&BD_success);
	if (!BD_success) {
		cout << "Error: PC_bsf_Init failed!" << endl;
		exit(1);
	};

	BD_success = true;
	BC_Init(&BD_success);
	if (!BD_success) {
		cout << "Error: BC_Init failed (not enough memory)!" << endl;
		exit(1);
	};

	BC_Master();

	return 0;
};

static void BC_Master() {// The head function of the master process.
	PC_bsf_ParametersOutput(BD_order.parameter);
	BD_iterCounter = 0;
	PC_bsf_SetInitParameter(&(BD_order.parameter));
	BD_t = -(double)time(NULL);
	do {
		PC_bsf_JobDispatcher(&(BD_order.parameter), &BD_newJobCase, &BD_exit);
		if (BD_exit) break;
		BD_jobCase = BD_newJobCase;
		if (BD_jobCase > PP_BSF_MAX_JOB_CASE) {
			cout << "BC_Master:Error: Job Case = " << BD_jobCase << " > PP_BSF_MAX_JOB_CASE = " << PP_BSF_MAX_JOB_CASE << endl;
			BD_exit = BD_EXIT;
			break;
		};
		BC_MasterMap(!BD_EXIT);
		BC_WorkerMap();
		BC_WorkerReduce();

		switch (BD_jobCase) {
			case 0:
			PC_bsf_ProcessResults(
				&BD_extendedReduceResult_P->elem, 
				BD_extendedReduceResult_P->reduceCounter,
				&(BD_order.parameter), 
				&BD_newJobCase, &BD_exit);
#ifdef PP_BSF_ITER_OUTPUT
			if (BD_iterCounter % PP_BSF_TRACE_COUNT == 0)
				PC_bsf_IterOutput(
					&BD_extendedReduceResult_P->elem, 
					BD_extendedReduceResult_P->reduceCounter, 
					BD_order.parameter,
					BD_t + (double)time(NULL),
					BD_newJobCase);
#endif // PP_BSF_ITER_OUTPUT
			break;
		case 1:
			PC_bsf_ProcessResults_1( 
				&BD_extendedReduceResult_P_1->elem, 
				BD_extendedReduceResult_P_1->reduceCounter,
				&(BD_order.parameter), 
				&BD_newJobCase, &BD_exit);
#ifdef PP_BSF_ITER_OUTPUT
			if (BD_iterCounter % PP_BSF_TRACE_COUNT == 0)
				PC_bsf_IterOutput_1(
					&BD_extendedReduceResult_P_1->elem, 
					BD_extendedReduceResult_P_1->reduceCounter, 
					BD_order.parameter,
					BD_t + (double)time(NULL),
					BD_newJobCase);
#endif // PP_BSF_ITER_OUTPUT
			break;
		case 2:
			PC_bsf_ProcessResults_2(
				&BD_extendedReduceResult_P_2->elem,
				BD_extendedReduceResult_P_2->reduceCounter,
				&(BD_order.parameter),
				&BD_newJobCase, &BD_exit);
#ifdef PP_BSF_ITER_OUTPUT
			if (BD_iterCounter % PP_BSF_TRACE_COUNT == 0)
				PC_bsf_IterOutput_2(
					&BD_extendedReduceResult_P_2->elem,
					BD_extendedReduceResult_P_2->reduceCounter,
					BD_order.parameter,
					BD_t + (double)time(NULL),
					BD_newJobCase);
#endif // PP_BSF_ITER_OUTPUT
			break;
		case 3:
			PC_bsf_ProcessResults_3(
				&BD_extendedReduceResult_P_3->elem,
				BD_extendedReduceResult_P_3->reduceCounter,
				&(BD_order.parameter),
				&BD_newJobCase, &BD_exit);
#ifdef PP_BSF_ITER_OUTPUT
			if (BD_iterCounter % PP_BSF_TRACE_COUNT == 0)
				PC_bsf_IterOutput_3(
					&BD_extendedReduceResult_P_3->elem,
					BD_extendedReduceResult_P_3->reduceCounter,
					BD_order.parameter,
					BD_t + (double)time(NULL),
					BD_newJobCase);
#endif // PP_BSF_ITER_OUTPUT
			break;
		default:
			cout << "BC_Master: Undefined job type!" << endl;
			break;
		}

		BD_iterCounter++;
	} while (!BD_exit);

	BD_t += (double)time(NULL);

	switch (BD_jobCase) {
	case 0:
		PC_bsf_ProblemOutput(&BD_extendedReduceResult_P->elem, BD_extendedReduceResult_P->reduceCounter, BD_order.parameter, BD_t);
		break;
	case 1:
		PC_bsf_ProblemOutput_1(&BD_extendedReduceResult_P_1->elem, BD_extendedReduceResult_P_1->reduceCounter, BD_order.parameter, BD_t);
		break;
	case 2:
		PC_bsf_ProblemOutput_2(&BD_extendedReduceResult_P_2->elem, BD_extendedReduceResult_P_2->reduceCounter, BD_order.parameter, BD_t);
		break;
	case 3:
		PC_bsf_ProblemOutput_3(&BD_extendedReduceResult_P_3->elem, BD_extendedReduceResult_P_3->reduceCounter, BD_order.parameter, BD_t);
		break;
	default:
		cout << "BC_Master: Undefined job type!" << endl;
		break;
	};
};

static void BC_MasterMap(bool exit) { // Forms an order and sends it to the worker processes to perform the Map function in the current iteration.
	PC_bsfAssignJobCase(BD_jobCase);
	PC_bsfAssignIterCounter(BD_iterCounter);

	BD_order.exit = exit;
	BD_order.jobCase = BD_jobCase;
	BD_order.iterCounter = BD_iterCounter;
}

static bool BC_WorkerMap() { // Performs the Map function
	if (BD_order.exit)
		return BD_EXIT;

	PC_bsfAssignJobCase(BD_order.jobCase);
	PC_bsfAssignIterCounter(BD_order.iterCounter);
	PC_bsfAssignSublistLength(BD_listSize);
	PC_bsfAssignAddressOffset(0);
	PC_bsfAssignParameter(BD_order.parameter);

#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
#pragma omp parallel for num_threads(PP_BSF_NUM_THREADS)
#else
#pragma omp parallel for
#endif // PP_BSF_NUM_THREADS
#endif // PP_BSF_OMP
	for (int i = 0; i < BD_listSize; i++) {

		PC_bsfAssignNumberInSublist(i);
		switch (BD_order.jobCase) {
		case 0:
			BD_extendedReduceList[i].reduceCounter = 1;
			PC_bsf_MapF(&BD_mapList[i], &BD_extendedReduceList[i].elem,
				&BD_extendedReduceList[i].reduceCounter);
			break;
		case 1:
			BD_extendedReduceList_1[i].reduceCounter = 1;
			PC_bsf_MapF_1(&BD_mapList[i], &BD_extendedReduceList_1[i].elem,
				&BD_extendedReduceList_1[i].reduceCounter);
			break;
		case 2:
			BD_extendedReduceList_2[i].reduceCounter = 1;
			PC_bsf_MapF_2(&BD_mapList[i], &BD_extendedReduceList_2[i].elem,
				&BD_extendedReduceList_2[i].reduceCounter);
			break;
		case 3:
			BD_extendedReduceList_3[i].reduceCounter = 1;
			PC_bsf_MapF_3(&BD_mapList[i], &BD_extendedReduceList_3[i].elem,
				&BD_extendedReduceList_3[i].reduceCounter);
			break;
		default:
			cout << "BC_WorkerReduce Error: Undefined job type!" << endl;
			break;
		};
	};
	return !BD_EXIT;
};

static void BC_WorkerReduce() {

	switch (BD_order.jobCase) {
	case 0:
		BC_ProcessExtendedReduceList(BD_extendedReduceList, BD_listSize, &BD_extendedReduceResult_P);
		break;
	case 1:
		BC_ProcessExtendedReduceList_1(BD_extendedReduceList_1, BD_listSize, &BD_extendedReduceResult_P_1);
		break;
	case 2:
		BC_ProcessExtendedReduceList_2(BD_extendedReduceList_2, BD_listSize, &BD_extendedReduceResult_P_2);
		break;
	case 3:
		BC_ProcessExtendedReduceList_3(BD_extendedReduceList_3, BD_listSize, &BD_extendedReduceResult_P_3);
		break;
	default:
		cout << "BC_WorkerReduce Error: Undefined job type!" << endl;
		break;
	};
};

static void BC_ProcessExtendedReduceList(BT_extendedReduceElem_T* reduceList, int length, BT_extendedReduceElem_T** extendedReduceResult_P) {
	int firstSuccessIndex = -1;

	*extendedReduceResult_P = &reduceList[0];

	for (int i = 0; i < length; i++) {
		if (reduceList[i].reduceCounter > 0) {
			*extendedReduceResult_P = &reduceList[i];
			firstSuccessIndex = i;
			break;
		};
	};

	if (firstSuccessIndex >= 0) {
		for (int i = firstSuccessIndex + 1; i < length; i++)
			if (BD_extendedReduceList[i].reduceCounter > 0) {
				PC_bsf_ReduceF(&(*extendedReduceResult_P)->elem, &BD_extendedReduceList[i].elem,
					&(*extendedReduceResult_P)->elem);
				(*extendedReduceResult_P)->reduceCounter += BD_extendedReduceList[i].reduceCounter;
			};
	};
};

static void BC_ProcessExtendedReduceList_1(BT_extendedReduceElem_T_1* reduceList, int length, BT_extendedReduceElem_T_1** extendedReduceResult_P) {
	int firstSuccessIndex = -1;

	*extendedReduceResult_P = &reduceList[0];

	for (int i = 0; i < length; i++) {
		if (reduceList[i].reduceCounter > 0) {
			*extendedReduceResult_P = &reduceList[i];
			firstSuccessIndex = i;
			break;
		};
	};

	if (firstSuccessIndex >= 0) {
		for (int i = firstSuccessIndex + 1; i < length; i++)
			if (BD_extendedReduceList_1[i].reduceCounter > 0) {
				PC_bsf_ReduceF_1(&(*extendedReduceResult_P)->elem, &BD_extendedReduceList_1[i].elem,
					&(*extendedReduceResult_P)->elem);
				(*extendedReduceResult_P)->reduceCounter += BD_extendedReduceList_1[i].reduceCounter;
			};
	};
};

static void BC_ProcessExtendedReduceList_2(BT_extendedReduceElem_T_2* reduceList, int length, BT_extendedReduceElem_T_2** extendedReduceResult_P) {
	int firstSuccessIndex = -1;

	*extendedReduceResult_P = &reduceList[0];

	for (int i = 0; i < length; i++) {
		if (reduceList[i].reduceCounter > 0) {
			*extendedReduceResult_P = &reduceList[i];
			firstSuccessIndex = i;
			break;
		};
	};

	if (firstSuccessIndex >= 0) {
		for (int i = firstSuccessIndex + 1; i < length; i++)
			if (BD_extendedReduceList_2[i].reduceCounter > 0) {
				PC_bsf_ReduceF_2(&(*extendedReduceResult_P)->elem, &BD_extendedReduceList_2[i].elem,
					&(*extendedReduceResult_P)->elem);
				(*extendedReduceResult_P)->reduceCounter += BD_extendedReduceList_2[i].reduceCounter;
			};
	};
};

static void BC_ProcessExtendedReduceList_3(BT_extendedReduceElem_T_3* reduceList, int length, BT_extendedReduceElem_T_3** extendedReduceResult_P) {
	int firstSuccessIndex = -1;

	*extendedReduceResult_P = &reduceList[0];

	for (int i = 0; i < length; i++) {
		if (reduceList[i].reduceCounter > 0) {
			*extendedReduceResult_P = &reduceList[i];
			firstSuccessIndex = i;
			break;
		};
	};

	if (firstSuccessIndex >= 0) {
		for (int i = firstSuccessIndex + 1; i < length; i++)
			if (BD_extendedReduceList_3[i].reduceCounter > 0) {
				PC_bsf_ReduceF_3(&(*extendedReduceResult_P)->elem, &BD_extendedReduceList_3[i].elem,
					&(*extendedReduceResult_P)->elem);
				(*extendedReduceResult_P)->reduceCounter += BD_extendedReduceList_3[i].reduceCounter;
			};
	};
};

static void BC_Init(bool* success) {// Performs the memory allocation and the initialization of the skeleton data structures and variables.

	BD_masterRank = 0;
	BD_numOfWorkers = 1;
	
	PC_bsf_SetListSize(&BD_listSize);

	PC_bsf_SetInitParameter(&(BD_order.parameter));

	for (int i = 0; i < BD_listSize; i++)
		PC_bsf_SetMapListElem(&BD_mapList[i], i);
}