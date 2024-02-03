/*==============================================================================
Project: LiFe
Theme: LPP Generator (no MPI)
Module: BSF-Forwards.h (Problem Independent Function Forwards)
Author: Leonid B. Sokolinsky 
This source code is a part of BSF Skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "BSF-Types.h"					// Problem Independent Types 

/*---------------- Assigning Values to BSF Problem Variables ----------------------*/
void PC_bsfAssignAddressOffset(int value);
void PC_bsfAssignIterCounter(int value);
void PC_bsfAssignJobCase(int value);
void PC_bsfAssignMpiMaster(int value);
void PC_bsfAssignMpiRank(int value);
void PC_bsfAssignNumberInSublist(int value);
void PC_bsfAssignNumOfWorkers(int value);
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter);
void PC_bsfAssignSublistLength(int value);

/*---------------- BSF Skeleton Functions ----------------------*/
static void	BC_Init(bool* success);		// Performs the memory allocation and the initialization of the skeleton data structures and variables.
static void BC_Master();				// The head function of the master process.
static void BC_MasterMap(bool exit);	// Forms an order and sends it to the worker processes to perform the Map function in the current iteration.
static void BC_ProcessExtendedReduceList(BT_extendedReduceElem_T* reduceList, int length, BT_extendedReduceElem_T** extendedReduceResult_P);
static void BC_ProcessExtendedReduceList_1(BT_extendedReduceElem_T_1* reduceList, int length, BT_extendedReduceElem_T_1** extendedReduceResult_P);
static void BC_ProcessExtendedReduceList_2(BT_extendedReduceElem_T_2* reduceList, int length, BT_extendedReduceElem_T_2** extendedReduceResult_P);
static void BC_ProcessExtendedReduceList_3(BT_extendedReduceElem_T_3* reduceList, int length, BT_extendedReduceElem_T_3** extendedReduceResult_P);
static bool BC_WorkerMap();				// Performs the Map function.
static void BC_WorkerReduce();			// Sends to the master process the element that is the sum of all reduce-sublist elements.