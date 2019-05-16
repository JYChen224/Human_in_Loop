/*
 * data1.h
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

#ifndef RTW_HEADER_data1_h_
#define RTW_HEADER_data1_h_
#include <string.h>
#ifndef data1_COMMON_INCLUDES_
# define data1_COMMON_INCLUDES_
#include <brtenv.h>
#include <rtkernel.h>
#include <rti_assert.h>
#include <rtidefineddatatypes.h>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* data1_COMMON_INCLUDES_ */

#include "data1_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rt_zcfcn.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->odeF = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

/* Block signals (auto storage) */
typedef struct {
  real_T SFunction1;                   /* '<S1>/S-Function1' */
  real_T Gain;                         /* '<Root>/Gain' */
  real_T freq;                         /* '<Root>/freq' */
  real_T Integrator;                   /* '<S9>/Integrator' */
  real_T Product4;                     /* '<S9>/Product4' */
  real_T Sum2;                         /* '<S9>/Sum2' */
  real_T In1;                          /* '<S4>/In1' */
  real_T real_weight;                  /* '<Root>/MATLAB Function' */
} B_data1_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T orig_weight;                  /* '<Root>/MATLAB Function' */
  boolean_T orig_weight_not_empty;     /* '<Root>/MATLAB Function' */
} DW_data1_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T Integrator_CSTATE;            /* '<S9>/Integrator' */
} X_data1_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T Integrator_CSTATE;            /* '<S9>/Integrator' */
} XDot_data1_T;

/* State disabled  */
typedef struct {
  boolean_T Integrator_CSTATE;         /* '<S9>/Integrator' */
} XDis_data1_T;

/* Zero-crossing (trigger) state */
typedef struct {
  ZCSigState TriggeredSubsystem_Trig_ZCE;/* '<Root>/Triggered Subsystem' */
} PrevZCX_data1_T;

#ifndef ODE1_INTG
#define ODE1_INTG

/* ODE1 Integration Data */
typedef struct {
  real_T *f[1];                        /* derivatives */
} ODE1_IntgData;

#endif

/* Parameters (auto storage) */
struct P_data1_T_ {
  real_T voltage_Y0;                   /* Computed Parameter: voltage_Y0
                                        * Referenced by: '<S4>/voltage'
                                        */
  real_T Button_Value;                 /* Expression: 0
                                        * Referenced by: '<Root>/Button'
                                        */
  real_T Gain_Gain;                    /* Expression: 10
                                        * Referenced by: '<Root>/Gain'
                                        */
  real_T current_weight_Value;         /* Expression: 0
                                        * Referenced by: '<Root>/current_weight'
                                        */
  real_T freq_Value;                   /* Expression: 1
                                        * Referenced by: '<Root>/freq'
                                        */
  real_T Integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S9>/Integrator'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_data1_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_data1_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeF[1][1];
  ODE1_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block parameters (auto storage) */
extern P_data1_T data1_P;

/* Block signals (auto storage) */
extern B_data1_T data1_B;

/* Continuous states (auto storage) */
extern X_data1_T data1_X;

/* Block states (auto storage) */
extern DW_data1_T data1_DW;

/* External data declarations for dependent source files */

/* Zero-crossing (trigger) state */
extern PrevZCX_data1_T data1_PrevZCX;

/* Model entry point functions */
extern void data1_initialize(void);
extern void data1_output(void);
extern void data1_update(void);
extern void data1_terminate(void);

/* Real-time Model object */
extern RT_MODEL_data1_T *const data1_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'data1'
 * '<S1>'   : 'data1/ADC_CLASS2_BL2'
 * '<S2>'   : 'data1/MATLAB Function'
 * '<S3>'   : 'data1/RTI Data'
 * '<S4>'   : 'data1/Triggered Subsystem'
 * '<S5>'   : 'data1/Varying Lowpass Filter'
 * '<S6>'   : 'data1/RTI Data/RTI Data Store'
 * '<S7>'   : 'data1/RTI Data/RTI Data Store/RTI Data Store'
 * '<S8>'   : 'data1/RTI Data/RTI Data Store/RTI Data Store/RTI Data Store'
 * '<S9>'   : 'data1/Varying Lowpass Filter/FOS'
 */
#endif                                 /* RTW_HEADER_data1_h_ */
