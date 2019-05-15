/***************************************************************************

   Source file Human_in_Loop_trc_ptr.c:

   Definition of function that initializes the global TRC pointers

   RTI1202 7.9 (02-Nov-2017)
   Sat May 11 07:01:28 2019

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
volatile int32_T *p_0_Human_in_Loop_int32_T_1 = NULL;
volatile uint32_T *p_0_Human_in_Loop_uint32_T_2 = NULL;
volatile uint8_T *p_0_Human_in_Loop_uint8_T_3 = NULL;
volatile boolean_T *p_0_Human_in_Loop_boolean_T_4 = NULL;
volatile real_T *p_0_Human_in_Loop_real_T_5 = NULL;
volatile real_T *p_0_Human_in_Loop_real_T_6 = NULL;
volatile real_T *p_0_Human_in_Loop_real_T_7 = NULL;
volatile real_T *p_0_Human_in_Loop_real_T_8 = NULL;
volatile real_T *p_0_Human_in_Loop_real_T_9 = NULL;
volatile real_T *p_0_Human_in_Loop_real_T_10 = NULL;
volatile real_T *p_0_Human_in_Loop_real_T_11 = NULL;
volatile real_T *p_0_Human_in_Loop_real_T_12 = NULL;
volatile real_T *p_0_Human_in_Loop_real_T_13 = NULL;
volatile real_T *p_0_Human_in_Loop_real_T_14 = NULL;
volatile real_T *p_1_Human_in_Loop_real_T_0 = NULL;
volatile int32_T *p_1_Human_in_Loop_int32_T_1 = NULL;
volatile int8_T *p_1_Human_in_Loop_int8_T_2 = NULL;
volatile boolean_T *p_1_Human_in_Loop_boolean_T_3 = NULL;
volatile real_T *p_1_Human_in_Loop_real_T_4 = NULL;
volatile uint32_T *p_1_Human_in_Loop_uint32_T_5 = NULL;
volatile boolean_T *p_1_Human_in_Loop_boolean_T_6 = NULL;
volatile real_T *p_2_Human_in_Loop_real_T_4 = NULL;
volatile uint32_T *p_2_Human_in_Loop_uint32_T_5 = NULL;
volatile int8_T *p_2_Human_in_Loop_int8_T_6 = NULL;
volatile boolean_T *p_2_Human_in_Loop_boolean_T_7 = NULL;
volatile real_T *p_2_Human_in_Loop_real_T_8 = NULL;
volatile real_T *p_2_Human_in_Loop_real_T_9 = NULL;
volatile real_T *p_2_Human_in_Loop_real_T_10 = NULL;
volatile real_T *p_2_Human_in_Loop_real_T_11 = NULL;
volatile real_T *p_2_Human_in_Loop_real_T_12 = NULL;
volatile real_T *p_2_Human_in_Loop_real_T_13 = NULL;
volatile real_T *p_3_Human_in_Loop_real_T_0 = NULL;

/*
 *  Declare the functions, that initially assign TRC pointers
 */
static void rti_init_trc_pointers_0(void);

/* Global pointers to data type transitions are separated in different functions to avoid overloading */
static void rti_init_trc_pointers_0(void)
{
  p_0_Human_in_Loop_real_T_0 = &Human_in_Loop_B.SFunction1;
  p_0_Human_in_Loop_int32_T_1 = &Human_in_Loop_B.SFunction1_o3;
  p_0_Human_in_Loop_uint32_T_2 = &Human_in_Loop_B.SFunction1_o2_b;
  p_0_Human_in_Loop_uint8_T_3 = &Human_in_Loop_B.SFunction1_o1_l[0];
  p_0_Human_in_Loop_boolean_T_4 = &Human_in_Loop_B.SFunction1_a;
  p_0_Human_in_Loop_real_T_5 = &Human_in_Loop_B.sf_Mux.x[0];
  p_0_Human_in_Loop_real_T_6 = &Human_in_Loop_B.sf_Mux_p.x[0];
  p_0_Human_in_Loop_real_T_7 = &Human_in_Loop_B.sf_MATLABFunction_l.y;
  p_0_Human_in_Loop_real_T_8 = &Human_in_Loop_B.sf_MATLABFunction_ns.y;
  p_0_Human_in_Loop_real_T_9 = &Human_in_Loop_B.sf_MATLABFunction_m.y;
  p_0_Human_in_Loop_real_T_10 = &Human_in_Loop_B.sf_MATLABFunction_g.y;
  p_0_Human_in_Loop_real_T_11 = &Human_in_Loop_B.sf_MATLABFunction_p.y;
  p_0_Human_in_Loop_real_T_12 = &Human_in_Loop_B.sf_MATLABFunction_g5.y;
  p_0_Human_in_Loop_real_T_13 = &Human_in_Loop_B.sf_Mux3.x[0];
  p_0_Human_in_Loop_real_T_14 = &Human_in_Loop_B.sf_Mux1_c.x[0];
  p_1_Human_in_Loop_real_T_0 = &Human_in_Loop_P.Controller_CALIB_SPEED;
  p_1_Human_in_Loop_int32_T_1 = &Human_in_Loop_P.LRN_time_delay;
  p_1_Human_in_Loop_int8_T_2 = &Human_in_Loop_P.OptTrigger_STRIDE_NUM;
  p_1_Human_in_Loop_boolean_T_3 = &Human_in_Loop_P.StateMachine_BT_CALIB;
  p_1_Human_in_Loop_real_T_4 = &Human_in_Loop_P.Mean_Y0;
  p_1_Human_in_Loop_uint32_T_5 = &Human_in_Loop_P.Delay1_DelayLength;
  p_1_Human_in_Loop_boolean_T_6 = &Human_in_Loop_P.VCC1_Value;
  p_2_Human_in_Loop_real_T_4 = &Human_in_Loop_DW.u4low2_states[0];
  p_2_Human_in_Loop_uint32_T_5 = &Human_in_Loop_DW.state[0];
  p_2_Human_in_Loop_int8_T_6 = &Human_in_Loop_DW.RT4_write_buf;
  p_2_Human_in_Loop_boolean_T_7 = &Human_in_Loop_DW.iter_not_empty;
  p_2_Human_in_Loop_real_T_8 = &Human_in_Loop_DW.sf_MATLABFunction_l.data_mem[0];
  p_2_Human_in_Loop_real_T_9 = &Human_in_Loop_DW.sf_MATLABFunction_ns.data_mem[0];
  p_2_Human_in_Loop_real_T_10 = &Human_in_Loop_DW.sf_MATLABFunction_m.data_mem[0];
  p_2_Human_in_Loop_real_T_11 = &Human_in_Loop_DW.sf_MATLABFunction_g.data_mem[0];
  p_2_Human_in_Loop_real_T_12 = &Human_in_Loop_DW.sf_MATLABFunction_p.data_mem[0];
  p_2_Human_in_Loop_real_T_13 = &Human_in_Loop_DW.sf_MATLABFunction_g5.data_mem
    [0];
  p_3_Human_in_Loop_real_T_0 = &Human_in_Loop_X.low_pass_CSTATE[0];
}

void Human_in_Loop_rti_init_trc_pointers(void)
{
  rti_init_trc_pointers_0();
}
