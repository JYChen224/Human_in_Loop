/*
 * Human_in_Loop_private.h
 *
 * Code generation for model "Human_in_Loop".
 *
 * Model version              : 1.987
 * Simulink Coder version : 8.13 (R2017b) 24-Jul-2017
 * C source code generated on : Wed May  8 14:27:10 2019
 *
 * Target selection: rti1202.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Custom Processor->Custom
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_Human_in_Loop_private_h_
#define RTW_HEADER_Human_in_Loop_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "Human_in_Loop.h"

extern DacCl1AnalogOutSDrvObject *pRTIDacC1AnalogOut_Ch_16;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_6;
extern DioCl2EncoderInSDrvObject *pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1;
extern DioCl2EncoderInSDrvObject *pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5;
extern DioCl1DigInSDrvObject *pRTIDioC1DigIn_Port_1_Ch_2;
extern DioCl1DigOutSDrvObject *pRTIDioC1DigOut_Port_3_Ch_11;
extern DioCl1DigOutSDrvObject *pRTIDioC1DigOut_Port_3_Ch_13;
extern DioCl1DigOutSDrvObject *pRTIDioC1DigOut_Port_1_Ch_1;
extern real_T rt_powd_snf(real_T u0, real_T u1);
extern real_T rt_hypotd_snf(real_T u0, real_T u1);
extern void Human_in_Loop_Mux(real_T rtu_x1, real_T rtu_x2,
  B_Mux_Human_in_Loop_T *localB);
extern void Human_in_Loo_ControlModule_Init(void);
extern void Human_in_Lo_ControlModule_Reset(void);
extern void Human_in_Loop_ControlModule(void);
extern void Human_in_Loo_ControlModule_Term(void);

#endif                                 /* RTW_HEADER_Human_in_Loop_private_h_ */
