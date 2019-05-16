/*********************** dSPACE target specific file *************************

   Header file data1_trc_ptr.h:

   Declaration of function that initializes the global TRC pointers

   RTI1202 7.9 (02-Nov-2017)
   Wed May 15 17:48:50 2019

   Copyright 2019, dSPACE GmbH. All rights reserved.

 *****************************************************************************/
#ifndef RTI_HEADER_data1_trc_ptr_h_
#define RTI_HEADER_data1_trc_ptr_h_

/* Include the model header file. */
#include "data1.h"
#include "data1_private.h"
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
EXTERN_C volatile real_T *p_0_data1_real_T_0;
EXTERN_C volatile real_T *p_1_data1_real_T_0;
EXTERN_C volatile real_T *p_2_data1_real_T_0;
EXTERN_C volatile boolean_T *p_2_data1_boolean_T_1;
EXTERN_C volatile real_T *p_3_data1_real_T_0;

/*
 *  Declare the general function for TRC pointer initialization
 */
EXTERN_C void data1_rti_init_trc_pointers(void);

#endif                                 /* RTI_HEADER_data1_trc_ptr_h_ */
