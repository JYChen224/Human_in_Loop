/*********************** dSPACE target specific file *************************

   Include file Human_in_Loop_rti.c:

   Definition of functions and variables for the system I/O and for
   the hardware and software interrupts used.

   RTI1202 7.9 (02-Nov-2017)
   Wed May  8 14:27:10 2019

   Copyright 2019, dSPACE GmbH. All rights reserved.

 *****************************************************************************/

#if !(defined(__RTI_SIMENGINE__) || defined(RTIMP_FRAME))
# error This file may be included only by the RTI(-MP) simulation engine.
#endif

/* Include the model header file. */
#include "Human_in_Loop.h"
#include "Human_in_Loop_private.h"

/* Defines for block output and parameter structure existence */
#define RTI_rtB_STRUCTURE_EXISTS       1
#define RTI_rtP_STRUCTURE_EXISTS       1
#define RTB_STRUCTURE_NAME             Human_in_Loop_B
#define RTP_STRUCTURE_NAME             Human_in_Loop_P

/* dSPACE generated includes for header files */
#include <brtenv.h>
#include <rtkernel.h>
#include <rti_assert.h>
#include <rtidefineddatatypes.h>
#ifndef dsRtmGetNumSampleTimes
# define dsRtmGetNumSampleTimes(rtm)   2
#endif

#ifndef dsRtmSetTaskTime
#define dsRtmSetTaskTime(rtm, sti, val) (((rtm)->Timing.taskTime0) = (val))
#endif

/****** Definitions: task functions for timer tasks *********************/

/* Timer Task 1. (Base rate). */
static void rti_TIMERA(rtk_p_task_control_block task)
{
  /* Task entry code BEGIN */
  /* -- None. -- */
  /* Task entry code END */

  /* Task code. */
  baseRateService(task);

  /* Task exit code BEGIN */
  /* -- None. -- */
  /* Task exit code END */
}

/****** Definitions: task functions for timer interrupts ****************/

/* Timer Interrupt: <Root>/Timer Interrupt */
static void rti_TIMERB(rtk_p_task_control_block task)
{
  /* Task entry code BEGIN */
  /* -- None. -- */
  /* Task entry code END */

  /* Task code. */
  {
    int32_T i;

    /* RateTransition: '<S4>/RT1' */
    switch (Human_in_Loop_DW.RT1_write_buf) {
     case 0:
      Human_in_Loop_DW.RT1_read_buf = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT1_read_buf = 0;
      break;

     default:
      Human_in_Loop_DW.RT1_read_buf = Human_in_Loop_DW.RT1_last_buf_wr;
      break;
    }

    if (Human_in_Loop_DW.RT1_read_buf != 0) {
      Human_in_Loop_B.RT1[0] = Human_in_Loop_DW.RT1_Buffer1[0];
      Human_in_Loop_B.RT1[1] = Human_in_Loop_DW.RT1_Buffer1[1];
    } else {
      Human_in_Loop_B.RT1[0] = Human_in_Loop_DW.RT1_Buffer0[0];
      Human_in_Loop_B.RT1[1] = Human_in_Loop_DW.RT1_Buffer0[1];
    }

    Human_in_Loop_DW.RT1_read_buf = -1;

    /* End of RateTransition: '<S4>/RT1' */

    /* RateTransition: '<S4>/RT2' */
    switch (Human_in_Loop_DW.RT2_write_buf) {
     case 0:
      Human_in_Loop_DW.RT2_read_buf = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT2_read_buf = 0;
      break;

     default:
      Human_in_Loop_DW.RT2_read_buf = Human_in_Loop_DW.RT2_last_buf_wr;
      break;
    }

    if (Human_in_Loop_DW.RT2_read_buf != 0) {
      Human_in_Loop_B.RT2[0] = Human_in_Loop_DW.RT2_Buffer1[0];
      Human_in_Loop_B.RT2[1] = Human_in_Loop_DW.RT2_Buffer1[1];
    } else {
      Human_in_Loop_B.RT2[0] = Human_in_Loop_DW.RT2_Buffer0[0];
      Human_in_Loop_B.RT2[1] = Human_in_Loop_DW.RT2_Buffer0[1];
    }

    Human_in_Loop_DW.RT2_read_buf = -1;

    /* End of RateTransition: '<S4>/RT2' */

    /* RateTransition: '<S4>/RT3' */
    switch (Human_in_Loop_DW.RT3_write_buf) {
     case 0:
      Human_in_Loop_DW.RT3_read_buf = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT3_read_buf = 0;
      break;

     default:
      Human_in_Loop_DW.RT3_read_buf = Human_in_Loop_DW.RT3_last_buf_wr;
      break;
    }

    if (Human_in_Loop_DW.RT3_read_buf != 0) {
      Human_in_Loop_B.RT3[0] = Human_in_Loop_DW.RT3_Buffer1[0];
      Human_in_Loop_B.RT3[1] = Human_in_Loop_DW.RT3_Buffer1[1];
      Human_in_Loop_B.RT3[2] = Human_in_Loop_DW.RT3_Buffer1[2];
    } else {
      Human_in_Loop_B.RT3[0] = Human_in_Loop_DW.RT3_Buffer0[0];
      Human_in_Loop_B.RT3[1] = Human_in_Loop_DW.RT3_Buffer0[1];
      Human_in_Loop_B.RT3[2] = Human_in_Loop_DW.RT3_Buffer0[2];
    }

    Human_in_Loop_DW.RT3_read_buf = -1;

    /* End of RateTransition: '<S4>/RT3' */

    /* RateTransition: '<S5>/RT1' */
    switch (Human_in_Loop_DW.RT1_write_buf_n) {
     case 0:
      Human_in_Loop_DW.RT1_read_buf_l = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT1_read_buf_l = 0;
      break;

     default:
      Human_in_Loop_DW.RT1_read_buf_l = Human_in_Loop_DW.RT1_last_buf_wr_m;
      break;
    }

    if (Human_in_Loop_DW.RT1_read_buf_l != 0) {
      Human_in_Loop_B.RT1_n[0] = Human_in_Loop_DW.RT1_Buffer1_f[0];
      Human_in_Loop_B.RT1_n[1] = Human_in_Loop_DW.RT1_Buffer1_f[1];
      Human_in_Loop_B.RT1_n[2] = Human_in_Loop_DW.RT1_Buffer1_f[2];
      Human_in_Loop_B.RT1_n[3] = Human_in_Loop_DW.RT1_Buffer1_f[3];
    } else {
      Human_in_Loop_B.RT1_n[0] = Human_in_Loop_DW.RT1_Buffer0_a[0];
      Human_in_Loop_B.RT1_n[1] = Human_in_Loop_DW.RT1_Buffer0_a[1];
      Human_in_Loop_B.RT1_n[2] = Human_in_Loop_DW.RT1_Buffer0_a[2];
      Human_in_Loop_B.RT1_n[3] = Human_in_Loop_DW.RT1_Buffer0_a[3];
    }

    Human_in_Loop_DW.RT1_read_buf_l = -1;

    /* End of RateTransition: '<S5>/RT1' */

    /* RateTransition: '<S2>/RT1' */
    switch (Human_in_Loop_DW.RT1_write_buf_e) {
     case 0:
      Human_in_Loop_DW.RT1_read_buf_p = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT1_read_buf_p = 0;
      break;

     default:
      Human_in_Loop_DW.RT1_read_buf_p = Human_in_Loop_DW.RT1_last_buf_wr_a;
      break;
    }

    for (i = 0; i < 5; i++) {
      if (Human_in_Loop_DW.RT1_read_buf_p != 0) {
        Human_in_Loop_B.RT1_a[i] = Human_in_Loop_DW.RT1_Buffer1_e[i];
      } else {
        Human_in_Loop_B.RT1_a[i] = Human_in_Loop_DW.RT1_Buffer0_o[i];
      }
    }

    Human_in_Loop_DW.RT1_read_buf_p = -1;

    /* End of RateTransition: '<S2>/RT1' */

    /* RateTransition: '<S2>/RT2' */
    switch (Human_in_Loop_DW.RT2_write_buf_i) {
     case 0:
      Human_in_Loop_DW.RT2_read_buf_k = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT2_read_buf_k = 0;
      break;

     default:
      Human_in_Loop_DW.RT2_read_buf_k = Human_in_Loop_DW.RT2_last_buf_wr_d;
      break;
    }

    if (Human_in_Loop_DW.RT2_read_buf_k != 0) {
      Human_in_Loop_B.RT2_h[0] = Human_in_Loop_DW.RT2_Buffer1_i[0];
      Human_in_Loop_B.RT2_h[1] = Human_in_Loop_DW.RT2_Buffer1_i[1];
      Human_in_Loop_B.RT2_h[2] = Human_in_Loop_DW.RT2_Buffer1_i[2];
      Human_in_Loop_B.RT2_h[3] = Human_in_Loop_DW.RT2_Buffer1_i[3];
    } else {
      Human_in_Loop_B.RT2_h[0] = Human_in_Loop_DW.RT2_Buffer0_k[0];
      Human_in_Loop_B.RT2_h[1] = Human_in_Loop_DW.RT2_Buffer0_k[1];
      Human_in_Loop_B.RT2_h[2] = Human_in_Loop_DW.RT2_Buffer0_k[2];
      Human_in_Loop_B.RT2_h[3] = Human_in_Loop_DW.RT2_Buffer0_k[3];
    }

    Human_in_Loop_DW.RT2_read_buf_k = -1;

    /* End of RateTransition: '<S2>/RT2' */

    /* S-Function (rti_commonblock): '<S6>/S-Function1' */
    Human_in_Loop_ControlModule();

    /* End of Outputs for S-Function (rti_commonblock): '<S6>/S-Function1' */
  }

  /* Task exit code BEGIN */
  /* -- None. -- */
  /* Task exit code END */
}

/* ===== Declarations of RTI blocks ======================================== */
DacCl1AnalogOutSDrvObject *pRTIDacC1AnalogOut_Ch_16;
AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_6;
DioCl2EncoderInSDrvObject *pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1;
DioCl2EncoderInSDrvObject *pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5;
DioCl1DigInSDrvObject *pRTIDioC1DigIn_Port_1_Ch_2;
DioCl1DigOutSDrvObject *pRTIDioC1DigOut_Port_3_Ch_11;
DioCl1DigOutSDrvObject *pRTIDioC1DigOut_Port_3_Ch_13;
DioCl1DigOutSDrvObject *pRTIDioC1DigOut_Port_1_Ch_1;

/* ===== Definition of interface functions for simulation engine =========== */
#if GRTINTERFACE == 1
#ifdef MULTITASKING
# define dsIsSampleHit(RTM,sti)        rtmGetSampleHitPtr(RTM)[sti]
#else
# define dsIsSampleHit(RTM,sti)        1
#endif

#else
#ifndef rtmStepTask
# define rtmStepTask(rtm, idx)         ((rtm)->Timing.TaskCounters.TID[(idx)] == 0)
#endif

# define dsIsSampleHit(RTM,sti)        rtmStepTask(RTM, sti)
#endif

#undef __INLINE
#if defined(_INLINE)
# define __INLINE                      static inline
#else
# define __INLINE                      static
#endif

/*Define additional variables*/
static time_T dsTFinal = -1.0;

#define dsGetTFinal(rtm)               (dsTFinal)

static time_T dsStepSize = 0.0002;

# define dsGetStepSize(rtm)            (dsStepSize)

static void rti_mdl_initialize_host_services(void)
{
  DsDaq_Init(0, 32, 1);
}

static void rti_mdl_initialize_io_boards(void)
{
  /* Registering of RTI products and modules at VCM */
  {
    vcm_module_register(VCM_MID_RTI1202, (void *) 0,
                        VCM_TXT_RTI1202, 7, 9, 0,
                        VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);

    {
      vcm_module_descriptor_type* msg_mod_ptr;
      msg_mod_ptr = vcm_module_register(VCM_MID_MATLAB, (void *) 0,
        VCM_TXT_MATLAB, 9, 3, 0,
        VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);
      vcm_module_register(VCM_MID_SIMULINK, msg_mod_ptr,
                          VCM_TXT_SIMULINK, 9, 0, 0,
                          VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);
      vcm_module_register(VCM_MID_RTW, msg_mod_ptr,
                          VCM_TXT_RTW, 8, 13, 0,
                          VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);
      vcm_module_register(VCM_MID_STATEFLOW, msg_mod_ptr,
                          VCM_TXT_STATEFLOW, 9, 0, 0,
                          VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);
      vcm_module_register(VCM_MID_STATEFLOW_CODER, msg_mod_ptr,
                          VCM_TXT_STATEFLOW_CODER, 8, 13, 0,
                          VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);
    }
  }

  /* --- Human_in_Loop/Control Module/Motor/DAC_CLASS1_BL1 --- */
  /* --- [RTI120X, DAC C1] - Channel: 16 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorCode = IOLIB_NO_ERROR;

    /* Init DAC CL1 AnalogOut driver object pRTIDacC1AnalogOut_Ch_16 */
    ioErrorCode = DacCl1AnalogOut_create(&pRTIDacC1AnalogOut_Ch_16,
      DAC_CLASS1_MASK_CH_16);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Apply all parameter settings. Finalize the construction of the AnalogOut driver object */
    ioErrorCode = DacCl1AnalogOut_apply(pRTIDacC1AnalogOut_Ch_16);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL1 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 1 - Port: 1 - Channel: 1 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorRetValue = IOLIB_NO_ERROR;

    /* Create EMC Encoder driver object pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1 */
    ioErrorRetValue = DioCl2EncoderIn_create
      (&pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1, DIO_CLASS2_CHANNEL_1,
       DIO_ENC_INSTANCE_1,DIO_ENC_IUSAGE_DISABLED);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* -- Position measurement enabled -- */

    /* Sets number of encoder lines for the selected incremental Encoder */
    ioErrorRetValue = DioCl2EncoderIn_setEncoderLines
      (pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1, (UInt32) 250000.0);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Sets minimum position for the selected incremental Encoder (in lines) */
    ioErrorRetValue = DioCl2EncoderIn_setMinPosition
      (pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1, (Float64) -125000);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Sets maximum position for the selected incremental Encoder (in lines) */
    ioErrorRetValue = DioCl2EncoderIn_setMaxPosition
      (pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1, (Float64) 125000);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Sets index signal usage mode for the selected incremental Encoder */
    ioErrorRetValue = DioCl2EncoderIn_setIndexMode
      (pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1, (Float64)
       DIO_ENC_IMODE_DISABLED);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* specify the minimum encoder velocity (in lines per second) */
    ioErrorRetValue = DioCl2EncoderIn_setMinEncVelocity
      (pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1, (Float64) 1);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Sets gated mode for the selected incremental Encoder */
    ioErrorRetValue = DioCl2EncoderIn_setGatedMode
      (pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1, DIO_ENC_IUSAGE_ENABLED);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorRetValue = DioCl2EncoderIn_setMeasurementInterval
      (pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1, 1);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Apply all parameter settings. Finalize the construction of the Encoder Obj driver object */
    ioErrorRetValue = DioCl2EncoderIn_apply
      (pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL3 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 3 - Port: 1 - Channel: 5 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorRetValue = IOLIB_NO_ERROR;

    /* Create EMC Encoder driver object pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5 */
    ioErrorRetValue = DioCl2EncoderIn_create
      (&pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5, DIO_CLASS2_CHANNEL_5,
       DIO_ENC_INSTANCE_3,DIO_ENC_IUSAGE_DISABLED);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* -- Position measurement enabled -- */

    /* Sets number of encoder lines for the selected incremental Encoder */
    ioErrorRetValue = DioCl2EncoderIn_setEncoderLines
      (pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5, (UInt32) 250000.0);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Sets minimum position for the selected incremental Encoder (in lines) */
    ioErrorRetValue = DioCl2EncoderIn_setMinPosition
      (pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5, (Float64) -125000);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Sets maximum position for the selected incremental Encoder (in lines) */
    ioErrorRetValue = DioCl2EncoderIn_setMaxPosition
      (pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5, (Float64) 125000);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Sets index signal usage mode for the selected incremental Encoder */
    ioErrorRetValue = DioCl2EncoderIn_setIndexMode
      (pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5, (Float64)
       DIO_ENC_IMODE_DISABLED);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* specify the minimum encoder velocity (in lines per second) */
    ioErrorRetValue = DioCl2EncoderIn_setMinEncVelocity
      (pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5, (Float64) 1);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Sets gated mode for the selected incremental Encoder */
    ioErrorRetValue = DioCl2EncoderIn_setGatedMode
      (pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5, DIO_ENC_IUSAGE_ENABLED);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorRetValue = DioCl2EncoderIn_setMeasurementInterval
      (pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5, 1);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Apply all parameter settings. Finalize the construction of the Encoder Obj driver object */
    ioErrorRetValue = DioCl2EncoderIn_apply
      (pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* --- Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_IN_BL1 --- */
  /* --- [RTI120X, BITIN] - Port: 1 - Channel: 2 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorCode = IOLIB_NO_ERROR;

    /* Init DIO CL1 DigIn driver object pRTIDioC1DigIn_Port_1_Ch_2 */
    ioErrorCode = DioCl1DigIn_create(&pRTIDioC1DigIn_Port_1_Ch_2,
      DIO_CLASS1_PORT_1, DIO_CLASS1_MASK_CH_2);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigIn_setInputFilter(pRTIDioC1DigIn_Port_1_Ch_2, 0);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigIn_apply(pRTIDioC1DigIn_Port_1_Ch_2);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigIn_start(pRTIDioC1DigIn_Port_1_Ch_2);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup1 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 11 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorCode = IOLIB_NO_ERROR;

    /* define variables required for BitOut initial functions */
    UInt32 outputDataInit = 0;

    /* Init DIO CL1 DigOut driver object pRTIDioC1DigOut_Port_3_Ch_11 */
    ioErrorCode = DioCl1DigOut_create(&pRTIDioC1DigOut_Port_3_Ch_11,
      DIO_CLASS1_PORT_3, DIO_CLASS1_MASK_CH_11);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigOut_setSignalVoltage(pRTIDioC1DigOut_Port_3_Ch_11,
      DIO_CLASS1_SIGNAL_5_0_V);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    outputDataInit = ( ( ( (UInt32)0) << (11 - 1)) | outputDataInit);

    /* write initialization value to digital output channel 11-11 on port 3 */
    ioErrorCode = DioCl1DigOut_setChMaskOutData(pRTIDioC1DigOut_Port_3_Ch_11,
      outputDataInit);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigOut_apply(pRTIDioC1DigOut_Port_3_Ch_11);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigOut_start(pRTIDioC1DigOut_Port_3_Ch_11);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup3 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 13 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorCode = IOLIB_NO_ERROR;

    /* define variables required for BitOut initial functions */
    UInt32 outputDataInit = 0;

    /* Init DIO CL1 DigOut driver object pRTIDioC1DigOut_Port_3_Ch_13 */
    ioErrorCode = DioCl1DigOut_create(&pRTIDioC1DigOut_Port_3_Ch_13,
      DIO_CLASS1_PORT_3, DIO_CLASS1_MASK_CH_13);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigOut_setSignalVoltage(pRTIDioC1DigOut_Port_3_Ch_13,
      DIO_CLASS1_SIGNAL_5_0_V);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    outputDataInit = ( ( ( (UInt32)0) << (13 - 1)) | outputDataInit);

    /* write initialization value to digital output channel 13-13 on port 3 */
    ioErrorCode = DioCl1DigOut_setChMaskOutData(pRTIDioC1DigOut_Port_3_Ch_13,
      outputDataInit);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigOut_apply(pRTIDioC1DigOut_Port_3_Ch_13);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigOut_start(pRTIDioC1DigOut_Port_3_Ch_13);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* --- Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_OUT_BL1 --- */
  /* --- [RTI120X, BITOUT] - Port: 1 - Channel: 1 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorCode = IOLIB_NO_ERROR;

    /* define variables required for BitOut initial functions */
    UInt32 outputDataInit = 0;

    /* Init DIO CL1 DigOut driver object pRTIDioC1DigOut_Port_1_Ch_1 */
    ioErrorCode = DioCl1DigOut_create(&pRTIDioC1DigOut_Port_1_Ch_1,
      DIO_CLASS1_PORT_1, DIO_CLASS1_MASK_CH_1);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigOut_setSignalVoltage(pRTIDioC1DigOut_Port_1_Ch_1,
      DIO_CLASS1_SIGNAL_5_0_V);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    outputDataInit = ( ( ( (UInt32)0) << (1 - 1)) | outputDataInit);

    /* write initialization value to digital output channel 1-1 on port 1 */
    ioErrorCode = DioCl1DigOut_setChMaskOutData(pRTIDioC1DigOut_Port_1_Ch_1,
      outputDataInit);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigOut_apply(pRTIDioC1DigOut_Port_1_Ch_1);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigOut_start(pRTIDioC1DigOut_Port_1_Ch_1);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* dSPACE I/O Board DS120XSTDADCC1 #0 */
  /* --- Human_in_Loop/Sensor Data/Torque module/ADC_CLASS1_BL6 --- */
  /* --- [RTI120X, ADC C1] - Channel: 6 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorCode = IOLIB_NO_ERROR;

    /* Init ADC CL1 AnalogIn driver object pRTIAdcC1AnalogIn_Ch_6 */
    ioErrorCode = AdcCl1AnalogIn_create(&pRTIAdcC1AnalogIn_Ch_6,
      ADC_CLASS1_CHANNEL_6);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Set parameters for ADC CL1 AnalogIn driver object pRTIAdcC1AnalogIn_Ch_6 */
    ioErrorCode = AdcCl1AnalogIn_setConversionMode(pRTIAdcC1AnalogIn_Ch_6,
      ADC_CLASS1_SINGLE_CONV_MODE);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = AdcCl1AnalogIn_setConversTrigger(pRTIAdcC1AnalogIn_Ch_6,
      ADC_CLASS1_TRIGGER_SW);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* dSPACE I/O Board DS120XSTDADCC1 #0 Unit:ADCC1 */

  /* dSPACE I/O Board DS120XSTDADCC1 #0 Unit:ADCC1 Group:ADC */
  /* --- Human_in_Loop/Sensor Data/Torque module/ADC_CLASS1_BL6 --- */
  /* --- [RTI120X, ADC C1] - Channel: 6 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorCode = IOLIB_NO_ERROR;

    /* Apply- and Start-Fcn for pRTIAdcC1AnalogIn_Ch_6 */
    ioErrorCode = AdcCl1AnalogIn_apply(pRTIAdcC1AnalogIn_Ch_6);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = AdcCl1AnalogIn_start(pRTIAdcC1AnalogIn_Ch_6);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }
}

/* Function rti_mdl_slave_load() is empty */
#define rti_mdl_slave_load()

/* Function rti_mdl_rtk_initialize() is empty */
#define rti_mdl_rtk_initialize()

static void rti_mdl_initialize_io_units(void)
{
  /* --- Human_in_Loop/Control Module/Motor/DAC_CLASS1_BL1 --- */
  /* --- [RTI120X, DAC C1] - Channel: 16 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorCode = IOLIB_NO_ERROR;

    /* All channel outputs are released from high impedance state */
    DacCl1AnalogOut_setOutputHighZ(pRTIDacC1AnalogOut_Ch_16,
      DAC_CLASS1_HIGH_Z_OFF);

    /* write initial value of CL1 DAC for output channel 16 */
    ioErrorCode = DacCl1AnalogOut_setOutputValue(pRTIDacC1AnalogOut_Ch_16,
      DAC_CLASS1_CHANNEL_16, (Float64) 0.0);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Activates AnalogOut functionality */
    ioErrorCode = DacCl1AnalogOut_start(pRTIDacC1AnalogOut_Ch_16);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL1 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 1 - Port: 1 - Channel: 1 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorRetValue = IOLIB_NO_ERROR;

    /* Resets start (initial) position for the selected incremental Encoder (in lines) */
    ioErrorRetValue = DioCl2EncoderIn_setEncPosition
      (pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1, (Float64) 0);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Activate Encoder signals read */
    ioErrorRetValue = DioCl2EncoderIn_start
      (pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL3 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 3 - Port: 1 - Channel: 5 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorRetValue = IOLIB_NO_ERROR;

    /* Resets start (initial) position for the selected incremental Encoder (in lines) */
    ioErrorRetValue = DioCl2EncoderIn_setEncPosition
      (pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5, (Float64) 0);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    /* Activate Encoder signals read */
    ioErrorRetValue = DioCl2EncoderIn_start
      (pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5);
    if (ioErrorRetValue > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup1 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 11 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorCode = IOLIB_NO_ERROR;

    /* define variables required for BitOut initial functions */
    UInt32 outputDataInit = 0;
    ioErrorCode = DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_3_Ch_11,
      DIO_CLASS1_HIGH_Z_OFF);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    outputDataInit = ( ( ( (UInt32)0) << (11 - 1)) | outputDataInit);

    /* write initialization value to digital output channel 11-11 on port 3 */
    ioErrorCode = DioCl1DigOut_setChMaskOutData(pRTIDioC1DigOut_Port_3_Ch_11,
      outputDataInit);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_11);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup3 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 13 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorCode = IOLIB_NO_ERROR;

    /* define variables required for BitOut initial functions */
    UInt32 outputDataInit = 0;
    ioErrorCode = DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_3_Ch_13,
      DIO_CLASS1_HIGH_Z_OFF);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    outputDataInit = ( ( ( (UInt32)0) << (13 - 1)) | outputDataInit);

    /* write initialization value to digital output channel 13-13 on port 3 */
    ioErrorCode = DioCl1DigOut_setChMaskOutData(pRTIDioC1DigOut_Port_3_Ch_13,
      outputDataInit);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_13);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }

  /* --- Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_OUT_BL1 --- */
  /* --- [RTI120X, BITOUT] - Port: 1 - Channel: 1 --- */
  {
    /* define a variable for IO error handling */
    UInt32 ioErrorCode = IOLIB_NO_ERROR;

    /* define variables required for BitOut initial functions */
    UInt32 outputDataInit = 0;
    ioErrorCode = DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_1_Ch_1,
      DIO_CLASS1_HIGH_Z_OFF);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    outputDataInit = ( ( ( (UInt32)0) << (1 - 1)) | outputDataInit);

    /* write initialization value to digital output channel 1-1 on port 1 */
    ioErrorCode = DioCl1DigOut_setChMaskOutData(pRTIDioC1DigOut_Port_1_Ch_1,
      outputDataInit);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }

    ioErrorCode = DioCl1DigOut_write(pRTIDioC1DigOut_Port_1_Ch_1);
    if (ioErrorCode > IOLIB_NO_ERROR) {
      RTLIB_EXIT(1);
    }
  }
}

/* Function rti_mdl_acknowledge_interrupts() is empty */
#define rti_mdl_acknowledge_interrupts()

/* Function rti_mdl_timetables_register() is empty */
#define rti_mdl_timetables_register()

/* Function rti_mdl_timesync_simstate() is empty */
#define rti_mdl_timesync_simstate()

/* Function rti_mdl_timesync_baserate() is empty */
#define rti_mdl_timesync_baserate()

static void rti_mdl_background(void)
{
  /* DsDaq background call */
  DsDaq_Background(0);
}

__INLINE void rti_mdl_sample_input(void)
{
  /* Calls for base sample time: [0.0002, 0.0] */
  /* dSPACE I/O Board DS120XSTDADCC1 #0 Unit:ADCC1 */

  /* dSPACE I/O Board DS120XSTDADCC1 #0 Unit:ADCC1 Group:ADC */
  /* fire burst- or conversion trigger for analog input channel. Called by ADC C1 block pRTIAdcC1AnalogIn_Ch_6 */
  AdcCl1AnalogIn_setConversSwTrigger(pRTIAdcC1AnalogIn_Ch_6);
  AdcCl1AnalogIn_write(pRTIAdcC1AnalogIn_Ch_6);

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL1 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 1 - Port: 1 - Channel: 1 --- */
  {
    /* define variables required for encoder sensor realtime functions */
    Float64 positionOrAngle= 0.0, speedOrAngVelocity= 0.0;
    UInt32 isIndexRaised= 0;

    /* Reads complete input data from selected encoder input channels (update data from hardware) */
    DioCl2EncoderIn_read(pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1);

    /* Gets encoder position (line) */
    DioCl2EncoderIn_getEncPosition(pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1,
      &positionOrAngle);

    /* Gets encoder speed (lines/second) */
    DioCl2EncoderIn_getEncVelocity(pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1,
      &speedOrAngVelocity);
    Human_in_Loop_B.SFunction1_o1 = (real_T) positionOrAngle;
    Human_in_Loop_B.SFunction1_o2 = (real_T) speedOrAngVelocity;
  }

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL3 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 3 - Port: 1 - Channel: 5 --- */
  {
    /* define variables required for encoder sensor realtime functions */
    Float64 positionOrAngle= 0.0, speedOrAngVelocity= 0.0;
    UInt32 isIndexRaised= 0;

    /* Reads complete input data from selected encoder input channels (update data from hardware) */
    DioCl2EncoderIn_read(pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5);

    /* Gets encoder position (line) */
    DioCl2EncoderIn_getEncPosition(pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5,
      &positionOrAngle);

    /* Gets encoder speed (lines/second) */
    DioCl2EncoderIn_getEncVelocity(pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5,
      &speedOrAngVelocity);
    Human_in_Loop_B.SFunction1_o1_k = (real_T) positionOrAngle;
    Human_in_Loop_B.SFunction1_o2_k = (real_T) speedOrAngVelocity;
  }

  /* --- Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_IN_BL1 --- */
  /* --- [RTI120X, BITIN] - Port: 1 - Channel: 2 --- */
  {
    UInt32 inputDataUInt32;

    /* get digital signal state on channel 2-2 on port 1 */
    DioCl1DigIn_read(pRTIDioC1DigIn_Port_1_Ch_2);
    DioCl1DigIn_getChMaskInData(pRTIDioC1DigIn_Port_1_Ch_2, &inputDataUInt32);
    Human_in_Loop_B.SFunction1_a = (boolean_T)((inputDataUInt32 &
      DIO_CLASS1_MASK_CH_2) >> 1);
  }
}

static void rti_mdl_daq_service()
{
  /* dSPACE Host Service */
  DsDaq_Service(0, 0, 1, (DsDaqSTimestampStruct *)
                rtk_current_task_absolute_time_ptr_get());
}

#undef __INLINE

/****** [EOF] ****************************************************************/
