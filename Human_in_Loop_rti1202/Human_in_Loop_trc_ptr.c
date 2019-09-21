/***************************************************************************

   Source file Human_in_Loop_trc_ptr.c:

   Definition of function that initializes the global TRC pointers

   RTI1202 7.9 (02-Nov-2017)
   Sat Sep 21 19:43:10 2019

   Copyright 2019, dSPACE GmbH. All rights reserved.

 *****************************************************************************/

/* Include header file. */
#include "Human_in_Loop_trc_ptr.h"
#include "Human_in_Loop.h"
#include "Human_in_Loop_private.h"

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
volatile real_T *p_0_Human_in_Loop_real_T_0 = NULL;
volatile boolean_T *p_0_Human_in_Loop_boolean_T_1 = NULL;
volatile real_T *p_0_Human_in_Loop_real_T_2 = NULL;
volatile real_T *p_0_Human_in_Loop_real_T_3 = NULL;
volatile real_T *p_1_Human_in_Loop_real_T_0 = NULL;
volatile int32_T *p_1_Human_in_Loop_int32_T_1 = NULL;
volatile int8_T *p_1_Human_in_Loop_int8_T_2 = NULL;
volatile boolean_T *p_1_Human_in_Loop_boolean_T_3 = NULL;
volatile real_T *p_1_Human_in_Loop_real_T_4 = NULL;
volatile uint32_T *p_1_Human_in_Loop_uint32_T_5 = NULL;
volatile boolean_T *p_1_Human_in_Loop_boolean_T_6 = NULL;
volatile real_T *p_2_Human_in_Loop_real_T_0 = NULL;
volatile int_T *p_2_Human_in_Loop_int_T_1 = NULL;
volatile int8_T *p_2_Human_in_Loop_int8_T_2 = NULL;
volatile boolean_T *p_2_Human_in_Loop_boolean_T_3 = NULL;
volatile real_T *p_3_Human_in_Loop_real_T_0 = NULL;

/*
 *  Declare the functions, that initially assign TRC pointers
 */
static void rti_init_trc_pointers_0(void);

/* Global pointers to data type transitions are separated in different functions to avoid overloading */
static void rti_init_trc_pointers_0(void)
{
  p_0_Human_in_Loop_real_T_0 = &Human_in_Loop_B.SFunction1;
  p_0_Human_in_Loop_boolean_T_1 = &Human_in_Loop_B.SFunction1_o3_b;
  p_0_Human_in_Loop_real_T_2 = &Human_in_Loop_B.sf_Mux.x[0];
  p_0_Human_in_Loop_real_T_3 = &Human_in_Loop_B.sf_Mux_p.x[0];
  p_1_Human_in_Loop_real_T_0 = &Human_in_Loop_P.Controller_FOLLOW_SLACK_ANGLE;
  p_1_Human_in_Loop_int32_T_1 = &Human_in_Loop_P.LRN_time_delay;
  p_1_Human_in_Loop_int8_T_2 = &Human_in_Loop_P.Controller_MODE;
  p_1_Human_in_Loop_boolean_T_3 = &Human_in_Loop_P.StateMachine_BT_CALIB;
  p_1_Human_in_Loop_real_T_4 = &Human_in_Loop_P.Mean_Y0;
  p_1_Human_in_Loop_uint32_T_5 = &Human_in_Loop_P.Delay3_DelayLength;
  p_1_Human_in_Loop_boolean_T_6 = &Human_in_Loop_P.VCC1_Value;
  p_2_Human_in_Loop_real_T_0 = &Human_in_Loop_DW.u4low2_states[0];
  p_2_Human_in_Loop_int_T_1 = &Human_in_Loop_DW.SFunction1_IWORK[0];
  p_2_Human_in_Loop_int8_T_2 = &Human_in_Loop_DW.RT4_write_buf;
  p_2_Human_in_Loop_boolean_T_3 = &Human_in_Loop_DW.last_footstate_not_empty;
  p_3_Human_in_Loop_real_T_0 = &Human_in_Loop_X.low_pass_CSTATE[0];
}

void Human_in_Loop_rti_init_trc_pointers(void)
{
  rti_init_trc_pointers_0();
}
