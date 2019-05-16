/***************************************************************************

   Source file data1_trc_ptr.c:

   Definition of function that initializes the global TRC pointers

   RTI1202 7.9 (02-Nov-2017)
   Wed May 15 17:48:50 2019

   Copyright 2019, dSPACE GmbH. All rights reserved.

 *****************************************************************************/

/* Include header file. */
#include "data1_trc_ptr.h"
#include "data1.h"
#include "data1_private.h"

/* Compiler options to turn off optimization. */
#if !defined(DS_OPTIMIZE_INIT_TRC_POINTERS)
#ifdef _MCCPPC

#pragma options -nOt -nOr -nOi -nOx

#endif

#ifdef __GNUC__

#pragma GCC optimize ("O0")

#endif

#ifdef _MSC_VER

#pragma optimize ("", off)

#endif
#endif

/* Definition of Global pointers to data type transitions (for TRC-file access) */
volatile real_T *p_0_data1_real_T_0 = NULL;
volatile real_T *p_1_data1_real_T_0 = NULL;
volatile real_T *p_2_data1_real_T_0 = NULL;
volatile boolean_T *p_2_data1_boolean_T_1 = NULL;
volatile real_T *p_3_data1_real_T_0 = NULL;

/*
 *  Declare the functions, that initially assign TRC pointers
 */
static void rti_init_trc_pointers_0(void);

/* Global pointers to data type transitions are separated in different functions to avoid overloading */
static void rti_init_trc_pointers_0(void)
{
  p_0_data1_real_T_0 = &data1_B.SFunction1;
  p_1_data1_real_T_0 = &data1_P.voltage_Y0;
  p_2_data1_real_T_0 = &data1_DW.orig_weight;
  p_2_data1_boolean_T_1 = &data1_DW.orig_weight_not_empty;
  p_3_data1_real_T_0 = &data1_X.Integrator_CSTATE;
}

void data1_rti_init_trc_pointers(void)
{
  rti_init_trc_pointers_0();
}
