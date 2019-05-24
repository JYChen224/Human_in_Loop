/*
 * Human_in_Loop_private.h
 *
 * Code generation for model "Human_in_Loop".
 *
 * Model version              : 1.1163
 * Simulink Coder version : 8.13 (R2017b) 24-Jul-2017
 * C source code generated on : Wed May 22 18:33:58 2019
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

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

/* ...  variable for information on a CAN channel */
extern can_tp1_canChannel* can_type1_channel_M1_C1;

/* ...  variable for information on a CAN channel */
extern can_tp1_canChannel* can_type1_channel_M1_C2;

/* ... definition of message variable for the RTICAN blocks */
#define CANTP1_M1_NUMMSG               3

extern can_tp1_canMsg* can_type1_msg_M1[CANTP1_M1_NUMMSG];

/* ... variable for taskqueue error checking                  */
extern Int32 rtican_type1_tq_error[CAN_TYPE1_NUM_MODULES]
  [CAN_TYPE1_NUM_TASKQUEUES];

/* Declaration of user indices (CAN_Type1_M1) */
#define CANTP1_M1_C1_RX_STD_0X1        0
#define RX_C1_STD_0X1                  0
#undef RX_C1_STD_0X1
#define CANTP1_M1_C1_RX_STD_0X64       1
#define RX_C1_STD_0X64                 1
#undef RX_C1_STD_0X64
#define CANTP1_M1_C2_TX_STD_0X64       2
#define TX_C2_STD_0X64                 2
#undef TX_C2_STD_0X64

/* predefine needed TX-definition code to support TX-Custom code */
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M1_C1_STD;
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M1_C1_XTD;

/* predefine pointer to CAN message object for STD-Msg */
extern can_tp1_canMsg* CANTP1_RX_SPMSG_M1_C1_STD;

/* predefine needed identifiers for mcr-queue */
extern can_tp1_canMsg* CANTP1_RX_M1_C1_MCRCLT_STD;

/* predefine needed TX-definition code to support TX-Custom code */
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M1_C2_STD;
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M1_C2_XTD;

/* predefine pointer to CAN message object for STD-Msg */
extern can_tp1_canMsg* CANTP1_RX_SPMSG_M1_C2_STD;

/* predefine needed identifiers for mcr-queue */
extern can_tp1_canMsg* CANTP1_RX_M1_C2_MCRCLT_STD;
extern DacCl1AnalogOutSDrvObject *pRTIDacC1AnalogOut_Ch_16;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_6;
extern DioCl2EncoderInSDrvObject *pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1;
extern DioCl2EncoderInSDrvObject *pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5;
extern DioCl1DigInSDrvObject *pRTIDioC1DigIn_Port_1_Ch_2;
extern DioCl1DigOutSDrvObject *pRTIDioC1DigOut_Port_3_Ch_11;
extern DioCl1DigOutSDrvObject *pRTIDioC1DigOut_Port_3_Ch_13;
extern DioCl1DigOutSDrvObject *pRTIDioC1DigOut_Port_1_Ch_1;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_1;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_2;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_3;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_4;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_9;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_10;

/* dSPACE I/O Board DS1202SER #1 Unit:GENSER Group:SETUP */
extern dsserChannel *rtiDS1202SER_B1_Ser[2];
extern real_T rt_hypotd_snf(real_T u0, real_T u1);
extern real_T rt_powd_snf(real_T u0, real_T u1);
extern real_T rt_roundd_snf(real_T u);
extern void Human_in_Loop_Mux1(real_T rtu_x1, real_T rtu_x2, real_T rtu_x3,
  real_T rtu_x4, real_T rtu_x5, real_T rtu_x6, B_Mux1_Human_in_Loop_T *localB);
extern void Human_in_Lo_MATLABFunction_Init(DW_MATLABFunction_Human_in_Lo_T
  *localDW);
extern void Human_in_Loop_MATLABFunction(real_T rtu_data,
  B_MATLABFunction_Human_in_Loo_T *localB, DW_MATLABFunction_Human_in_Lo_T
  *localDW);
extern void Human_in_Loop_Mux(real_T rtu_x1, real_T rtu_x2,
  B_Mux_Human_in_Loop_T *localB);
extern void Human_SerialDecodingSystem_Init(void);
extern void Huma_SerialDecodingSystem_Reset(void);
extern void Human_in_L_SerialDecodingSystem(void);
extern void Human_in_Loo_ControlModule_Init(void);
extern void Human_in_Lo_ControlModule_Reset(void);
extern void Human_in_Loop_ControlModule(void);
extern void Human_in_Loo_ControlModule_Term(void);

/* private model entry point functions */
extern void Human_in_Loop_derivatives(void);

#endif                                 /* RTW_HEADER_Human_in_Loop_private_h_ */
