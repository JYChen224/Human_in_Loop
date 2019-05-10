/*
 * Human_in_Loop_types.h
 *
 * Code generation for model "Human_in_Loop".
 *
 * Model version              : 1.1206
 * Simulink Coder version : 8.13 (R2017b) 24-Jul-2017
 * C source code generated on : Sat May 11 05:28:40 2019
 *
 * Target selection: rti1202.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Custom Processor->Custom
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_Human_in_Loop_types_h_
#define RTW_HEADER_Human_in_Loop_types_h_
#include "rtwtypes.h"
#include "multiword_types.h"

/* Custom Type definition for MATLAB Function: '<S25>/Bayesian Function' */
#ifndef typedef_coder_internal_ref_Human_in_L_T
#define typedef_coder_internal_ref_Human_in_L_T

typedef struct {
  real_T contents;
} coder_internal_ref_Human_in_L_T;

#endif                                 /*typedef_coder_internal_ref_Human_in_L_T*/

#ifndef typedef_cell_1_Human_in_Loop_T
#define typedef_cell_1_Human_in_Loop_T

typedef struct {
  coder_internal_ref_Human_in_L_T *f1;
  real_T f2[60];
  real_T f3;
} cell_1_Human_in_Loop_T;

#endif                                 /*typedef_cell_1_Human_in_Loop_T*/

#ifndef typedef_coder_internal_nested_functio_T
#define typedef_coder_internal_nested_functio_T

typedef struct {
  cell_1_Human_in_Loop_T environment;
} coder_internal_nested_functio_T;

#endif                                 /*typedef_coder_internal_nested_functio_T*/

#ifndef typedef_cell_2_Human_in_Loop_T
#define typedef_cell_2_Human_in_Loop_T

typedef struct {
  coder_internal_nested_functio_T f1;
  real_T f2[900];
  real_T f3;
  real_T f4[30];
} cell_2_Human_in_Loop_T;

#endif                                 /*typedef_cell_2_Human_in_Loop_T*/

#ifndef typedef_cell_3_Human_in_Loop_T
#define typedef_cell_3_Human_in_Loop_T

typedef struct {
  coder_internal_ref_Human_in_L_T *f1;
  coder_internal_nested_functio_T f2;
  real_T f3[900];
  real_T f4;
} cell_3_Human_in_Loop_T;

#endif                                 /*typedef_cell_3_Human_in_Loop_T*/

/* Custom Type definition for MATLAB Function: '<S4>/CMAES' */
#ifndef struct_emxArray_real_T_1x6
#define struct_emxArray_real_T_1x6

struct emxArray_real_T_1x6
{
  real_T data[6];
  int32_T size[2];
};

#endif                                 /*struct_emxArray_real_T_1x6*/

#ifndef typedef_emxArray_real_T_1x6_Human_in__T
#define typedef_emxArray_real_T_1x6_Human_in__T

typedef struct emxArray_real_T_1x6 emxArray_real_T_1x6_Human_in__T;

#endif                                 /*typedef_emxArray_real_T_1x6_Human_in__T*/

#ifndef struct_emxArray_real_T_2x6
#define struct_emxArray_real_T_2x6

struct emxArray_real_T_2x6
{
  real_T data[12];
  int32_T size[2];
};

#endif                                 /*struct_emxArray_real_T_2x6*/

#ifndef typedef_emxArray_real_T_2x6_Human_in__T
#define typedef_emxArray_real_T_2x6_Human_in__T

typedef struct emxArray_real_T_2x6 emxArray_real_T_2x6_Human_in__T;

#endif                                 /*typedef_emxArray_real_T_2x6_Human_in__T*/

#ifndef struct_emxArray_real_T_2x2
#define struct_emxArray_real_T_2x2

struct emxArray_real_T_2x2
{
  real_T data[4];
  int32_T size[2];
};

#endif                                 /*struct_emxArray_real_T_2x2*/

#ifndef typedef_emxArray_real_T_2x2_Human_in__T
#define typedef_emxArray_real_T_2x2_Human_in__T

typedef struct emxArray_real_T_2x2 emxArray_real_T_2x2_Human_in__T;

#endif                                 /*typedef_emxArray_real_T_2x2_Human_in__T*/

#ifndef struct_emxArray_real_T_2
#define struct_emxArray_real_T_2

struct emxArray_real_T_2
{
  real_T data[2];
  int32_T size;
};

#endif                                 /*struct_emxArray_real_T_2*/

#ifndef typedef_emxArray_real_T_2_Human_in_Lo_T
#define typedef_emxArray_real_T_2_Human_in_Lo_T

typedef struct emxArray_real_T_2 emxArray_real_T_2_Human_in_Lo_T;

#endif                                 /*typedef_emxArray_real_T_2_Human_in_Lo_T*/

/* Parameters (auto storage) */
typedef struct P_Human_in_Loop_T_ P_Human_in_Loop_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_Human_in_Loop_T RT_MODEL_Human_in_Loop_T;

#endif                                 /* RTW_HEADER_Human_in_Loop_types_h_ */
