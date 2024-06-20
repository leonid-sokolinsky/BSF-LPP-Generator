/*==============================================================================
Project: LiFe
Theme: LPP Generator (MPI)
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
	char emptystring[] = "";
	char* message = emptystring;
	unsigned send_success;
	unsigned recv_success;
	bool success;
	BC_MpiRun();
	BD_success = true;
	PC_bsf_Init(&BD_success);
	send_success = (unsigned)BD_success;
	MPI_Allreduce(&send_success, &recv_success, 1, MPI_UNSIGNED, MPI_PROD, MPI_COMM_WORLD);
	success = (bool)recv_success;
	if (!success) {
		if (BD_rank == BD_masterRank)
			cout << "Error: PC_bsf_Init failed!" << endl;
		MPI_Finalize();
		exit(1);
	}
	BD_success = true;
	BC_Init(&BD_success);
	MPI_Allreduce(&BD_success, &success, 1, MPI_UNSIGNED, MPI_LAND, MPI_COMM_WORLD);
	if (!success) {
		if (BD_rank == BD_masterRank)
			cout << "Error: BC_Init failed (not enough memory)!" << endl;
		MPI_Finalize();
		exit(1);
	}

	if (BD_rank == BD_masterRank) {
		BC_Master();
	}
	else {
		BC_Worker();
	}
	MPI_Finalize();
	return 0;
}

static void BC_Master() {// The head function of the master process.
	PC_bsf_ParametersOutput(BD_order.parameter);
	BD_iterCounter = 0;
	PC_bsf_SetInitParameter(&(BD_order.parameter));
	BD_t = -MPI_Wtime();
	do {
		PC_bsf_JobDispatcher(&(BD_order.parameter), &BD_newJobCase, &BD_exit, BD_t + MPI_Wtime());
		if (BD_exit) break;
		BD_jobCase = BD_newJobCase;
		if (BD_jobCase > PP_BSF_MAX_JOB_CASE) {
			cout << "BC_Master:Error: Job Case = " << BD_jobCase << " > PP_BSF_MAX_JOB_CASE = " << PP_BSF_MAX_JOB_CASE << endl;
			BD_exit = BD_EXIT;
			break;
		}
		BC_MasterMap(!BD_EXIT);
		BC_MasterReduce();
		//
		//
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
					BD_t + MPI_Wtime(),
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
					BD_t + MPI_Wtime(),
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
					BD_t + MPI_Wtime(),
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
					BD_t + MPI_Wtime(),
					BD_newJobCase);
#endif // PP_BSF_ITER_OUTPUT
			break;
		default:
			cout << "BC_Master: Undefined job type!" << endl;
			break;
		}

		BD_iterCounter++;
	} while (!BD_exit);

	BD_t += MPI_Wtime();

	BC_MasterMap(BD_EXIT);

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
	}
}

static void BC_Worker() {// The head function of a worker process.
	bool exit;
	while (true) {
		exit = BC_WorkerMap();
		if (exit) break;
		BC_WorkerReduce();
	}
}

static void BC_MasterMap(bool exit) { // Forms an order and sends it to the worker processes to perform the Map function in the current iteration.
	PC_bsfAssignJobCase(BD_jobCase);
	PC_bsfAssignIterCounter(BD_iterCounter);

	BD_order.exit = exit;
	BD_order.jobCase = BD_jobCase;
	BD_order.iterCounter = BD_iterCounter;

	for (int rank = 0; rank < BD_numOfWorkers; rank++) {
		MPI_Isend(
			&BD_order,
			sizeof(BT_order_T),
			MPI_BYTE,
			rank,
			0,
			MPI_COMM_WORLD,
			&BD_request[rank]);
	}
	MPI_Waitall(BD_numOfWorkers, BD_request, BD_status);
}

static void BC_MasterReduce() {

	for (int rank = 0; rank < BD_numOfWorkers; rank++)
		switch (BD_jobCase) {
		case 0:
			MPI_Irecv(&BD_extendedReduceList[rank], sizeof(BT_extendedReduceElem_T), MPI_BYTE, rank, 0, MPI_COMM_WORLD, &BD_request[rank]);
			break;
		case 1:
			MPI_Irecv(&BD_extendedReduceList_1[rank], sizeof(BT_extendedReduceElem_T_1), MPI_BYTE, rank, 0, MPI_COMM_WORLD, &BD_request[rank]);
			break;
		case 2:
			MPI_Irecv(&BD_extendedReduceList_2[rank], sizeof(BT_extendedReduceElem_T_2), MPI_BYTE, rank, 0, MPI_COMM_WORLD, &BD_request[rank]);
			break;
		case 3:
			MPI_Irecv(&BD_extendedReduceList_3[rank], sizeof(BT_extendedReduceElem_T_3), MPI_BYTE, rank, 0, MPI_COMM_WORLD, &BD_request[rank]);
			break;
		default:
			cout << "BC_MasterReduce Error: Undefined job type!" << endl;
			break;
		}
	MPI_Waitall(BD_numOfWorkers, BD_request, BD_status);

	switch (BD_jobCase) {
	case 0:
		BC_ProcessExtendedReduceList(BD_extendedReduceList, BD_numOfWorkers, &BD_extendedReduceResult_P);
		break;
	case 1:
		BC_ProcessExtendedReduceList_1(BD_extendedReduceList_1, BD_numOfWorkers, &BD_extendedReduceResult_P_1);
		break;
	case 2:
		BC_ProcessExtendedReduceList_2(BD_extendedReduceList_2, BD_numOfWorkers, &BD_extendedReduceResult_P_2);
		break;
	case 3:
		BC_ProcessExtendedReduceList_3(BD_extendedReduceList_3, BD_numOfWorkers, &BD_extendedReduceResult_P_3);
		break;
	default:
		cout << "BC_MasterReduce Error: Undefined job type!" << endl;
		break;
	}
}

static bool BC_WorkerMap() { // Performs the Map function
	MPI_Recv(
		&BD_order,
		sizeof(BT_order_T),
		MPI_BYTE,
		BD_masterRank,
		0,
		MPI_COMM_WORLD,
		&BD_status[BD_rank]);

	if (BD_order.exit)
		return BD_EXIT;

	PC_bsfAssignJobCase(BD_order.jobCase);
	PC_bsfAssignIterCounter(BD_order.iterCounter);
	PC_bsfAssignSublistLength(BD_sublistSize[BD_rank]);
	PC_bsfAssignAddressOffset(BD_offset[BD_rank]);
	PC_bsfAssignParameter(BD_order.parameter);
	PC_bsf_MapInit(BD_order.parameter);
#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
#pragma omp parallel for num_threads(PP_BSF_NUM_THREADS)
#else
#pragma omp parallel for
#endif // PP_BSF_NUM_THREADS
#endif // PP_BSF_OMP
	for (int i = BD_offset[BD_rank]; i < BD_offset[BD_rank] + BD_sublistSize[BD_rank]; i++) {
		int index;

		index = i - BD_offset[BD_rank];
		PC_bsfAssignNumberInSublist(index);
		switch (BD_order.jobCase) {
		case 0:
			BD_extendedReduceList[index].reduceCounter = 1;
			PC_bsf_MapF(&BD_mapSubList[index], &BD_extendedReduceList[index].elem,
				&BD_extendedReduceList[index].reduceCounter);
			break;
		case 1:
			BD_extendedReduceList_1[index].reduceCounter = 1;
			PC_bsf_MapF_1(&BD_mapSubList[index], &BD_extendedReduceList_1[index].elem,
				&BD_extendedReduceList_1[index].reduceCounter);
			break;
		case 2:
			BD_extendedReduceList_2[index].reduceCounter = 1;
			PC_bsf_MapF_2(&BD_mapSubList[index], &BD_extendedReduceList_2[index].elem,
				&BD_extendedReduceList_2[index].reduceCounter);
			break;
		case 3:
			BD_extendedReduceList_3[index].reduceCounter = 1;
			PC_bsf_MapF_3(&BD_mapSubList[index], &BD_extendedReduceList_3[index].elem,
				&BD_extendedReduceList_3[index].reduceCounter);
			break;
		default:
			cout << "BC_WorkerReduce Error: Undefined job type!" << endl;
			break;
		}
	}
	return !BD_EXIT;
}

static void BC_WorkerReduce() {

	switch (BD_order.jobCase) {
	case 0:
		BC_ProcessExtendedReduceList(BD_extendedReduceList, BD_sublistSize[BD_rank], &BD_extendedReduceResult_P);
		MPI_Send(BD_extendedReduceResult_P, sizeof(BT_extendedReduceElem_T), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD);
		break;
	case 1:
		BC_ProcessExtendedReduceList_1(BD_extendedReduceList_1, BD_sublistSize[BD_rank], &BD_extendedReduceResult_P_1);
		MPI_Send(BD_extendedReduceResult_P_1, sizeof(BT_extendedReduceElem_T_1), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD);
		break;
	case 2:
		BC_ProcessExtendedReduceList_2(BD_extendedReduceList_2, BD_sublistSize[BD_rank], &BD_extendedReduceResult_P_2);
		MPI_Send(BD_extendedReduceResult_P_2, sizeof(BT_extendedReduceElem_T_2), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD);
		break;
	case 3:
		BC_ProcessExtendedReduceList_3(BD_extendedReduceList_3, BD_sublistSize[BD_rank], &BD_extendedReduceResult_P_3);
		MPI_Send(BD_extendedReduceResult_P_3, sizeof(BT_extendedReduceElem_T_3), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD);
		break;
	default:
		cout << "BC_WorkerReduce Error: Undefined job type!" << endl;
		break;
	}
}

static void BC_ProcessExtendedReduceList(BT_extendedReduceElem_T* reduceList, int length, BT_extendedReduceElem_T** extendedReduceResult_P) {
	int firstSuccessIndex = -1;

	*extendedReduceResult_P = &reduceList[0];

	for (int i = 0; i < length; i++) {
		if (reduceList[i].reduceCounter > 0) {
			*extendedReduceResult_P = &reduceList[i];
			firstSuccessIndex = i;
			break;
		}
	}

	if (firstSuccessIndex >= 0) {
		for (int i = firstSuccessIndex + 1; i < length; i++)
			if (BD_extendedReduceList[i].reduceCounter > 0) {
				PC_bsf_ReduceF(&(*extendedReduceResult_P)->elem, &BD_extendedReduceList[i].elem,
					&(*extendedReduceResult_P)->elem);
				(*extendedReduceResult_P)->reduceCounter += BD_extendedReduceList[i].reduceCounter;
			}
	}
}

static void BC_ProcessExtendedReduceList_1(BT_extendedReduceElem_T_1* reduceList, int length, BT_extendedReduceElem_T_1** extendedReduceResult_P) {
	int firstSuccessIndex = -1;

	*extendedReduceResult_P = &reduceList[0];

	for (int i = 0; i < length; i++) {
		if (reduceList[i].reduceCounter > 0) {
			*extendedReduceResult_P = &reduceList[i];
			firstSuccessIndex = i;
			break;
		}
	}

	if (firstSuccessIndex >= 0) {
		for (int i = firstSuccessIndex + 1; i < length; i++)
			if (BD_extendedReduceList_1[i].reduceCounter > 0) {
				PC_bsf_ReduceF_1(&(*extendedReduceResult_P)->elem, &BD_extendedReduceList_1[i].elem,
					&(*extendedReduceResult_P)->elem);
				(*extendedReduceResult_P)->reduceCounter += BD_extendedReduceList_1[i].reduceCounter;
			}
	}
}

static void BC_ProcessExtendedReduceList_2(BT_extendedReduceElem_T_2* reduceList, int length, BT_extendedReduceElem_T_2** extendedReduceResult_P) {
	int firstSuccessIndex = -1;

	*extendedReduceResult_P = &reduceList[0];

	for (int i = 0; i < length; i++) {
		if (reduceList[i].reduceCounter > 0) {
			*extendedReduceResult_P = &reduceList[i];
			firstSuccessIndex = i;
			break;
		}
	}

	if (firstSuccessIndex >= 0) {
		for (int i = firstSuccessIndex + 1; i < length; i++)
			if (BD_extendedReduceList_2[i].reduceCounter > 0) {
				PC_bsf_ReduceF_2(&(*extendedReduceResult_P)->elem, &BD_extendedReduceList_2[i].elem,
					&(*extendedReduceResult_P)->elem);
				(*extendedReduceResult_P)->reduceCounter += BD_extendedReduceList_2[i].reduceCounter;
			}
	}
}

static void BC_ProcessExtendedReduceList_3(BT_extendedReduceElem_T_3* reduceList, int length, BT_extendedReduceElem_T_3** extendedReduceResult_P) {
	int firstSuccessIndex = -1;

	*extendedReduceResult_P = &reduceList[0];

	for (int i = 0; i < length; i++) {
		if (reduceList[i].reduceCounter > 0) {
			*extendedReduceResult_P = &reduceList[i];
			firstSuccessIndex = i;
			break;
		}
	}

	if (firstSuccessIndex >= 0) {
		for (int i = firstSuccessIndex + 1; i < length; i++)
			if (BD_extendedReduceList_3[i].reduceCounter > 0) {
				PC_bsf_ReduceF_3(&(*extendedReduceResult_P)->elem, &BD_extendedReduceList_3[i].elem,
					&(*extendedReduceResult_P)->elem);
				(*extendedReduceResult_P)->reduceCounter += BD_extendedReduceList_3[i].reduceCounter;
			}
	}
}

static void BC_Init(bool* success) {// Performs the memory allocation and the initialization of the skeleton data structures and variables.
	cout << setprecision(PP_BSF_PRECISION);
	int offset, first, last, subListSize;
	///
	PC_bsf_SetListSize(&BD_listSize);
	if (BD_size > BD_listSize + 1) {
		if (BD_rank == 0) cout << "Error: MPI_SIZE must be < Map List Size + 2 =" << BD_listSize + 2 << endl;
		MPI_Finalize();
		exit(1);
	}
	BD_masterRank = BD_size - 1;
	BD_numOfWorkers = BD_size - 1;
	PC_bsf_SetInitParameter(&(BD_order.parameter));

	if (BD_rank == BD_masterRank) {
		subListSize = BD_numOfWorkers;
	}
	else {

		BD_elemsPerWorker = BD_listSize / BD_numOfWorkers;
		BD_tailLength = BD_listSize - BD_elemsPerWorker * BD_numOfWorkers;

		offset = 0;
		for (int rank = 0; rank < BD_numOfWorkers; rank++) {
			BD_sublistSize[rank] = (int)(BD_elemsPerWorker + (rank < BD_tailLength ? 1 : 0));
			BD_offset[rank] = offset;
			offset += BD_sublistSize[rank];
		}

		subListSize = BD_sublistSize[BD_rank];

		BD_mapSubList = (PT_bsf_mapElem_T*)calloc(subListSize, sizeof(PT_bsf_mapElem_T));
		if (BD_mapSubList == NULL) {
			*success = false;
			return;
		}

		first = BD_offset[BD_rank];
		last = BD_offset[BD_rank] + BD_sublistSize[BD_rank] - 1;
		for (int i = first; i <= last; i++)
			PC_bsf_SetMapListElem(&BD_mapSubList[i - first], i);
	}

	BD_extendedReduceList = (BT_extendedReduceElem_T*)calloc(subListSize, sizeof(BT_extendedReduceElem_T));
	if (BD_extendedReduceList == NULL) { *success = false; return; }

	if (PP_BSF_MAX_JOB_CASE > 0) {
		BD_extendedReduceList_1 = (BT_extendedReduceElem_T_1*)calloc(subListSize, sizeof(BT_extendedReduceElem_T_1));
		if (BD_extendedReduceList_1 == NULL) { *success = false; return; }
	}

	if (PP_BSF_MAX_JOB_CASE > 1) {
		BD_extendedReduceList_2 = (BT_extendedReduceElem_T_2*)calloc(subListSize, sizeof(BT_extendedReduceElem_T_2));
		if (BD_extendedReduceList_2 == NULL) { *success = false; return; }
	}

	if (PP_BSF_MAX_JOB_CASE > 2) {
		BD_extendedReduceList_3 = (BT_extendedReduceElem_T_3*)calloc(subListSize, sizeof(BT_extendedReduceElem_T_3));
		if (BD_extendedReduceList_3 == NULL) { *success = false; return; }
	}
}

static void BC_MpiRun() {
	int rc;
	rc = MPI_Init(NULL, NULL);	// Starting MPI
	if (rc != MPI_SUCCESS) {
		printf("Error starting MPI program. Terminating! \n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &BD_rank);
	PC_bsfAssignMpiRank(BD_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &BD_size);
	PC_bsfAssignMpiMaster(BD_size - 1);
	if (BD_size > PP_BSF_MAX_MPI_SIZE) {
		if (BD_rank == 0) cout << "Error: MPI_SIZE exceeded the maximum allowed value PP_BSF_MAX_MPI_SIZE = " << PP_BSF_MAX_MPI_SIZE << endl;
		MPI_Finalize();
		exit(1);
	}
	PC_bsfAssignNumOfWorkers(BD_size - 1);
	if (BD_size < 2) {
		if (BD_rank == 0) cout << "Error: MPI_SIZE must be > 1" << endl;
		MPI_Finalize();
		exit(1);
	}
}