/*********************** dSPACE target specific file *************************

   Header file Human_in_Loop_trc_ptr.h:

   Declaration of function that initializes the global TRC pointers

   RTI1202 7.9 (02-Nov-2017)
   Wed May  8 17:05:26 2019

   Copyright 2019, dSPACE GmbH. All rights reserved.

 *****************************************************************************/
#ifndef RTI_HEADER_Human_in_Loop_trc_ptr_h_
#define RTI_HEADER_Human_in_Loop_trc_ptr_h_

/* Include the model header file. */
#include "Human_in_Loop.h"
#include "Human_in_Loop_private.h"
#ifdef EXTERN_C
#undef EXTERN_C
#endif

#ifdef __cplusplus
#define EXTERN_C                       extern "C"
#else
#define EXTERN_C                       extern
#endif

/*
 *  Declare the global TRC pointers
 */
EXTERN_C volatile real_T *p_0_Human_in_Loop_real_T_0;
EXTERN_C volatile boolean_T *p_0_Human_in_Loop_boolean_T_1;
EXTERN_C volatile real_T *p_0_Human_in_Loop_real_T_2;
EXTERN_C volatile real_T *p_0_Human_in_Loop_real_T_3;
EXTERN_C volatile real_T *p_1_Human_in_Loop_real_T_0;
EXTERN_C volatile int32_T *p_1_Human_in_Loop_int32_T_1;
EXTERN_C volatile int8_T *p_1_Human_in_Loop_int8_T_2;
EXTERN_C volatile boolean_T *p_1_Human_in_Loop_boolean_T_3;
EXTERN_C volatile real_T *p_1_Human_in_Loop_real_T_4;
EXTERN_C volatile boolean_T *p_1_Human_in_Loop_boolean_T_5;
EXTERN_C volatile real_T *p_2_Human_in_Loop_real_T_0;
EXTERN_C volatile int8_T *p_2_Human_in_Loop_int8_T_1;

/*
 *  Declare the general function for TRC pointer initialization
 */
EXTERN_C void Human_in_Loop_rti_init_trc_pointers(void);

#endif                                 /* RTI_HEADER_Human_in_Loop_trc_ptr_h_ */
