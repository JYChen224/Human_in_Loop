/*
 * data1.c
 *
 * Code generation for model "data1".
 *
 * Model version              : 1.4
 * Simulink Coder version : 8.13 (R2017b) 24-Jul-2017
 * C source code generated on : Wed May 15 17:48:50 2019
 *
 * Target selection: rti1202.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Custom Processor->Custom
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "data1_trc_ptr.h"
#include "data1.h"
#include "data1_private.h"

/* Block signals (auto storage) */
B_data1_T data1_B;

/* Continuous states */
X_data1_T data1_X;

/* Block states (auto storage) */
DW_data1_T data1_DW;

/* Previous zero-crossings (trigger) states */
PrevZCX_data1_T data1_PrevZCX;

/* Real-time model */
RT_MODEL_data1_T data1_M_;
RT_MODEL_data1_T *const data1_M = &data1_M_;

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 1;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  data1_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void data1_output(void)
{
  ZCEventType zcEvent;
  if (rtmIsMajorTimeStep(data1_M)) {
    /* set solver stop time */
    if (!(data1_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&data1_M->solverInfo, ((data1_M->Timing.clockTickH0
        + 1) * data1_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&data1_M->solverInfo, ((data1_M->Timing.clockTick0 +
        1) * data1_M->Timing.stepSize0 + data1_M->Timing.clockTickH0 *
        data1_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(data1_M)) {
    data1_M->Timing.t[0] = rtsiGetT(&data1_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(data1_M)) {
    /* S-Function (rti_commonblock): '<S1>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* Gain: '<Root>/Gain' */
    data1_B.Gain = data1_P.Gain_Gain * data1_B.SFunction1;

    /* MATLAB Function: '<Root>/MATLAB Function' incorporates:
     *  Constant: '<Root>/Button'
     *  Constant: '<Root>/current_weight'
     */
    /* MATLAB Function 'MATLAB Function': '<S2>:1' */
    if (!data1_DW.orig_weight_not_empty) {
      /* '<S2>:1:3' */
      /* '<S2>:1:4' */
      data1_DW.orig_weight = data1_P.current_weight_Value;
      data1_DW.orig_weight_not_empty = true;
    }

    if (data1_P.Button_Value == 1.0) {
      /* '<S2>:1:6' */
      /* '<S2>:1:7' */
      data1_B.real_weight = data1_P.current_weight_Value;

      /* '<S2>:1:8' */
      data1_DW.orig_weight = data1_P.current_weight_Value;
    } else {
      /* '<S2>:1:10' */
      data1_B.real_weight = data1_DW.orig_weight;
    }

    /* End of MATLAB Function: '<Root>/MATLAB Function' */

    /* Constant: '<Root>/freq' */
    data1_B.freq = data1_P.freq_Value;
  }

  /* Integrator: '<S9>/Integrator' */
  data1_B.Integrator = data1_X.Integrator_CSTATE;

  /* Product: '<S9>/Product4' */
  data1_B.Product4 = data1_B.freq * data1_B.Integrator;

  /* Outputs for Triggered SubSystem: '<Root>/Triggered Subsystem' incorporates:
   *  TriggerPort: '<S4>/Trigger'
   */
  if (rtmIsMajorTimeStep(data1_M) && rtmIsMajorTimeStep(data1_M)) {
    /* Constant: '<Root>/Button' */
    zcEvent = rt_ZCFcn(RISING_ZERO_CROSSING,
                       &data1_PrevZCX.TriggeredSubsystem_Trig_ZCE,
                       (data1_P.Button_Value));
    if (zcEvent != NO_ZCEVENT) {
      /* Inport: '<S4>/In1' */
      data1_B.In1 = data1_B.Product4;
    }
  }

  /* End of Outputs for SubSystem: '<Root>/Triggered Subsystem' */

  /* Sum: '<S9>/Sum2' */
  data1_B.Sum2 = data1_B.Gain - data1_B.Product4;
}

/* Model update function */
void data1_update(void)
{
  if (rtmIsMajorTimeStep(data1_M)) {
    rt_ertODEUpdateContinuousStates(&data1_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++data1_M->Timing.clockTick0)) {
    ++data1_M->Timing.clockTickH0;
  }

  data1_M->Timing.t[0] = rtsiGetSolverStopTime(&data1_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The resolution of this integer timer is 0.001, which is the step size
     * of the task. Size of "clockTick1" ensures timer will not overflow during the
     * application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    data1_M->Timing.clockTick1++;
    if (!data1_M->Timing.clockTick1) {
      data1_M->Timing.clockTickH1++;
    }
  }
}

/* Derivatives for root system: '<Root>' */
void data1_derivatives(void)
{
  XDot_data1_T *_rtXdot;
  _rtXdot = ((XDot_data1_T *) data1_M->derivs);

  /* Derivatives for Integrator: '<S9>/Integrator' */
  _rtXdot->Integrator_CSTATE = data1_B.Sum2;
}

/* Model initialize function */
void data1_initialize(void)
{
  /* Registration code */

  /* initialize real-time model */
  (void) memset((void *)data1_M, 0,
                sizeof(RT_MODEL_data1_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&data1_M->solverInfo, &data1_M->Timing.simTimeStep);
    rtsiSetTPtr(&data1_M->solverInfo, &rtmGetTPtr(data1_M));
    rtsiSetStepSizePtr(&data1_M->solverInfo, &data1_M->Timing.stepSize0);
    rtsiSetdXPtr(&data1_M->solverInfo, &data1_M->derivs);
    rtsiSetContStatesPtr(&data1_M->solverInfo, (real_T **) &data1_M->contStates);
    rtsiSetNumContStatesPtr(&data1_M->solverInfo, &data1_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&data1_M->solverInfo,
      &data1_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&data1_M->solverInfo,
      &data1_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&data1_M->solverInfo,
      &data1_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&data1_M->solverInfo, (&rtmGetErrorStatus(data1_M)));
    rtsiSetRTModelPtr(&data1_M->solverInfo, data1_M);
  }

  rtsiSetSimTimeStep(&data1_M->solverInfo, MAJOR_TIME_STEP);
  data1_M->intgData.f[0] = data1_M->odeF[0];
  data1_M->contStates = ((X_data1_T *) &data1_X);
  rtsiSetSolverData(&data1_M->solverInfo, (void *)&data1_M->intgData);
  rtsiSetSolverName(&data1_M->solverInfo,"ode1");
  rtmSetTPtr(data1_M, &data1_M->Timing.tArray[0]);
  data1_M->Timing.stepSize0 = 0.001;

  /* block I/O */
  (void) memset(((void *) &data1_B), 0,
                sizeof(B_data1_T));

  /* states (continuous) */
  {
    (void) memset((void *)&data1_X, 0,
                  sizeof(X_data1_T));
  }

  /* states (dwork) */
  (void) memset((void *)&data1_DW, 0,
                sizeof(DW_data1_T));

  {
    /* user code (registration function declaration) */
    /*Initialize global TRC pointers. */
    data1_rti_init_trc_pointers();
  }

  data1_PrevZCX.TriggeredSubsystem_Trig_ZCE = UNINITIALIZED_ZCSIG;

  /* InitializeConditions for Integrator: '<S9>/Integrator' */
  data1_X.Integrator_CSTATE = data1_P.Integrator_IC;

  /* SystemInitialize for MATLAB Function: '<Root>/MATLAB Function' */
  data1_DW.orig_weight_not_empty = false;

  /* SystemInitialize for Triggered SubSystem: '<Root>/Triggered Subsystem' */
  /* SystemInitialize for Outport: '<S4>/voltage' */
  data1_B.In1 = data1_P.voltage_Y0;

  /* End of SystemInitialize for SubSystem: '<Root>/Triggered Subsystem' */
}

/* Model terminate function */
void data1_terminate(void)
{
  /* (no terminate code required) */
}
