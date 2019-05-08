/*
 * Human_in_Loop.h
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

#ifndef RTW_HEADER_Human_in_Loop_h_
#define RTW_HEADER_Human_in_Loop_h_
#include <math.h>
#include <string.h>
#ifndef Human_in_Loop_COMMON_INCLUDES_
# define Human_in_Loop_COMMON_INCLUDES_
#include <brtenv.h>
#include <rtkernel.h>
#include <rti_assert.h>
#include <rtidefineddatatypes.h>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* Human_in_Loop_COMMON_INCLUDES_ */

#include "Human_in_Loop_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rt_nonfinite.h"
#include "rtGetInf.h"

/* Macros for accessing real-time model data structure */
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
# define rtmGetT(rtm)                  ((rtm)->Timing.taskTime0)
#endif

/* Block signals for system '<S19>/Mux' */
typedef struct {
  real_T x[2];                         /* '<S19>/Mux' */
} B_Mux_Human_in_Loop_T;

/* Block signals (auto storage) */
typedef struct {
  real_T result_data[3000];
  real_T tmp_data[3000];
  real_T b_result_data[2250];
  real_T tmp_data_m[2250];
  real_T torque_track[750];
  real_T torque_delta_track[750];
  real_T tmp_data_c[750];
  real_T tmp_data_k[750];
  real_T tmp_data_cx[750];
  real_T z1_data[750];
  real_T z1_data_b[750];
  real_T SFunction1;                   /* '<S38>/S-Function1' */
  real_T Gain;                         /* '<S21>/Gain' */
  real_T u4low2;                       /* '<S21>/0.4low2' */
  real_T x2k1;                         /* '<S37>/Unit Delay1' */
  real_T x1k1;                         /* '<S37>/Unit Delay' */
  real_T Gain1;                        /* '<S37>/Gain1' */
  real_T Gain2;                        /* '<S37>/Gain2' */
  real_T UnitDelay2;                   /* '<S37>/Unit Delay2' */
  real_T Gain4;                        /* '<S37>/Gain4' */
  real_T Add2;                         /* '<S37>/Add2' */
  real_T Gain3;                        /* '<S37>/Gain3' */
  real_T x2k;                          /* '<S37>/Add1' */
  real_T RT1[2];                       /* '<S4>/RT1' */
  real_T SFunction1_o1;                /* '<S27>/S-Function1' */
  real_T SFunction1_o2;                /* '<S27>/S-Function1' */
  real_T Gain_l;                       /* '<S19>/Gain' */
  real_T Gain1_b;                      /* '<S19>/Gain1' */
  real_T RT2[2];                       /* '<S4>/RT2' */
  real_T SFunction1_o1_k;              /* '<S28>/S-Function1' */
  real_T SFunction1_o2_k;              /* '<S28>/S-Function1' */
  real_T Gain2_h;                      /* '<S19>/Gain2' */
  real_T Gain3_m;                      /* '<S19>/Gain3' */
  real_T x2k1_k;                       /* '<S23>/Unit Delay1' */
  real_T x1k1_m;                       /* '<S23>/Unit Delay' */
  real_T Gain1_m;                      /* '<S23>/Gain1' */
  real_T Gain2_c;                      /* '<S23>/Gain2' */
  real_T UnitDelay2_b;                 /* '<S23>/Unit Delay2' */
  real_T Gain4_g;                      /* '<S23>/Gain4' */
  real_T Add2_m;                       /* '<S23>/Add2' */
  real_T Gain3_k;                      /* '<S23>/Gain3' */
  real_T x2k_i;                        /* '<S23>/Add1' */
  real_T RT3[3];                       /* '<S4>/RT3' */
  real_T RT1_n[4];                     /* '<S5>/RT1' */
  real_T RT1_a[5];                     /* '<S2>/RT1' */
  real_T RT2_h[4];                     /* '<S2>/RT2' */
  real_T Gain_g;                       /* '<S37>/Gain' */
  real_T x1k;                          /* '<S37>/Add' */
  real_T TSamp;                        /* '<S26>/TSamp' */
  real_T Uk1;                          /* '<S26>/UD' */
  real_T Diff;                         /* '<S26>/Diff' */
  real_T u4low1;                       /* '<S19>/0.4low1' */
  real_T UnitDelay;                    /* '<S22>/Unit Delay' */
  real_T UnitDelay1;                   /* '<S22>/Unit Delay1' */
  real_T Add1;                         /* '<S22>/Add1' */
  real_T Gain_p;                       /* '<S22>/Gain' */
  real_T Gain1_d;                      /* '<S22>/Gain1' */
  real_T Add2_i;                       /* '<S22>/Add2' */
  real_T Add3;                         /* '<S22>/Add3' */
  real_T Gain2_i;                      /* '<S22>/Gain2' */
  real_T Gain_i;                       /* '<S23>/Gain' */
  real_T x1k_i;                        /* '<S23>/Add' */
  real_T Gain2_j;                      /* '<S8>/Gain2' */
  real_T Gain1_p;                      /* '<S8>/Gain1' */
  real_T torque_des;                   /* '<S1>/Torque track' */
  real_T torque_delta_des;             /* '<S1>/Torque track' */
  real_T torque_trace[1500];           /* '<S1>/Torque track' */
  real_T torque_delta_trace[1500];     /* '<S1>/Torque track' */
  real_T vel;                          /* '<S8>/MATLAB Function' */
  real_T motor_vel_cmd;                /* '<S1>/Controller' */
  real_T mode;                         /* '<S5>/State Machine' */
  real_T state;                        /* '<S5>/State Machine' */
  real_T stride_time;                  /* '<S5>/State Machine' */
  real_T stride_timer;                 /* '<S5>/State Machine' */
  real_T x[4];                         /* '<S5>/Mux1' */
  real_T torque_dot;                   /* '<S21>/MATLAB Function' */
  real_T torque;                       /* '<S21>/Data process' */
  real_T state_c;                      /* '<S20>/FootSwitch Filter' */
  real_T x_d[3];                       /* '<S19>/Mux2' */
  real_T x_i[3];                       /* '<S19>/Mux1' */
  real_T angle;                        /* '<S19>/Data process1' */
  real_T angle_m;                      /* '<S19>/Data process' */
  real_T x_k[4];                       /* '<S13>/Mux1' */
  real_T x_o[5];                       /* '<S12>/Mux1' */
  boolean_T SFunction1_a;              /* '<S34>/S-Function1' */
  B_Mux_Human_in_Loop_T sf_Mux;        /* '<S21>/Mux' */
  B_Mux_Human_in_Loop_T sf_Mux_p;      /* '<S19>/Mux' */
} B_Human_in_Loop_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T u4low2_states[3];             /* '<S21>/0.4low2' */
  real_T UnitDelay1_DSTATE;            /* '<S37>/Unit Delay1' */
  real_T UnitDelay_DSTATE;             /* '<S37>/Unit Delay' */
  real_T UnitDelay2_DSTATE;            /* '<S37>/Unit Delay2' */
  real_T UnitDelay1_DSTATE_h;          /* '<S23>/Unit Delay1' */
  real_T UnitDelay_DSTATE_e;           /* '<S23>/Unit Delay' */
  real_T UnitDelay2_DSTATE_a;          /* '<S23>/Unit Delay2' */
  real_T UD_DSTATE;                    /* '<S26>/UD' */
  real_T u4low1_states[3];             /* '<S19>/0.4low1' */
  real_T UnitDelay_DSTATE_j;           /* '<S22>/Unit Delay' */
  real_T UnitDelay1_DSTATE_a;          /* '<S22>/Unit Delay1' */
  real_T u4low2_tmp;                   /* '<S21>/0.4low2' */
  volatile real_T RT1_Buffer0[2];      /* '<S4>/RT1' */
  volatile real_T RT1_Buffer1[2];      /* '<S4>/RT1' */
  volatile real_T RT2_Buffer0[2];      /* '<S4>/RT2' */
  volatile real_T RT2_Buffer1[2];      /* '<S4>/RT2' */
  volatile real_T RT3_Buffer0[3];      /* '<S4>/RT3' */
  volatile real_T RT3_Buffer1[3];      /* '<S4>/RT3' */
  volatile real_T RT1_Buffer0_a[4];    /* '<S5>/RT1' */
  volatile real_T RT1_Buffer1_f[4];    /* '<S5>/RT1' */
  volatile real_T RT1_Buffer0_o[5];    /* '<S2>/RT1' */
  volatile real_T RT1_Buffer1_e[5];    /* '<S2>/RT1' */
  volatile real_T RT2_Buffer0_k[4];    /* '<S2>/RT2' */
  volatile real_T RT2_Buffer1_i[4];    /* '<S2>/RT2' */
  real_T u4low1_tmp;                   /* '<S19>/0.4low1' */
  real_T TorqueMem[4400];              /* '<Root>/Data Store Memory' */
  real_T last_footstate;               /* '<S1>/Torque track' */
  real_T calib_state;                  /* '<S1>/Controller' */
  real_T reg_stride_time;              /* '<S5>/State Machine' */
  real_T reg_stride_time_count;        /* '<S5>/State Machine' */
  real_T reg_mode;                     /* '<S5>/State Machine' */
  real_T reg_state;                    /* '<S5>/State Machine' */
  real_T bt_run;                       /* '<S5>/State Machine' */
  real_T reg_last_switch;              /* '<S5>/State Machine' */
  real_T data[15];                     /* '<S21>/MATLAB Function' */
  real_T torque_zero;                  /* '<S21>/Data process' */
  real_T foot_state;                   /* '<S20>/FootSwitch Filter' */
  real_T filter_time;                  /* '<S20>/FootSwitch Filter' */
  real_T angle_zero;                   /* '<S19>/Data process1' */
  real_T angle_zero_f;                 /* '<S19>/Data process' */
  volatile int8_T RT1_write_buf;       /* '<S4>/RT1' */
  volatile int8_T RT1_read_buf;        /* '<S4>/RT1' */
  volatile int8_T RT1_last_buf_wr;     /* '<S4>/RT1' */
  volatile int8_T RT2_write_buf;       /* '<S4>/RT2' */
  volatile int8_T RT2_read_buf;        /* '<S4>/RT2' */
  volatile int8_T RT2_last_buf_wr;     /* '<S4>/RT2' */
  volatile int8_T RT3_write_buf;       /* '<S4>/RT3' */
  volatile int8_T RT3_read_buf;        /* '<S4>/RT3' */
  volatile int8_T RT3_last_buf_wr;     /* '<S4>/RT3' */
  volatile int8_T RT1_write_buf_n;     /* '<S5>/RT1' */
  volatile int8_T RT1_read_buf_l;      /* '<S5>/RT1' */
  volatile int8_T RT1_last_buf_wr_m;   /* '<S5>/RT1' */
  volatile int8_T RT1_write_buf_e;     /* '<S2>/RT1' */
  volatile int8_T RT1_read_buf_p;      /* '<S2>/RT1' */
  volatile int8_T RT1_last_buf_wr_a;   /* '<S2>/RT1' */
  volatile int8_T RT2_write_buf_i;     /* '<S2>/RT2' */
  volatile int8_T RT2_read_buf_k;      /* '<S2>/RT2' */
  volatile int8_T RT2_last_buf_wr_d;   /* '<S2>/RT2' */
} DW_Human_in_Loop_T;

/* Parameters (auto storage) */
struct P_Human_in_Loop_T_ {
  real_T Controller_CALIB_SPEED;       /* Mask Parameter: Controller_CALIB_SPEED
                                        * Referenced by: '<S1>/Controller'
                                        */
  real_T Controller_CALIB_TORQUE;      /* Mask Parameter: Controller_CALIB_TORQUE
                                        * Referenced by: '<S1>/Controller'
                                        */
  real_T Controller_FOLLOW_SLACK_ANGLE;/* Mask Parameter: Controller_FOLLOW_SLACK_ANGLE
                                        * Referenced by: '<S1>/Controller'
                                        */
  real_T DiscreteDerivative_ICPrevScaled;/* Mask Parameter: DiscreteDerivative_ICPrevScaled
                                          * Referenced by: '<S26>/UD'
                                          */
  real_T MATLABFunction_MAX_MOTOR_ANGLE;/* Mask Parameter: MATLABFunction_MAX_MOTOR_ANGLE
                                         * Referenced by: '<S8>/MATLAB Function'
                                         */
  real_T MATLABFunction_MAX_SPEED;     /* Mask Parameter: MATLABFunction_MAX_SPEED
                                        * Referenced by: '<S8>/MATLAB Function'
                                        */
  real_T MATLABFunction_MAX_TORQUE;    /* Mask Parameter: MATLABFunction_MAX_TORQUE
                                        * Referenced by: '<S8>/MATLAB Function'
                                        */
  real_T MATLABFunction_MIN_MOTOR_ANGLE;/* Mask Parameter: MATLABFunction_MIN_MOTOR_ANGLE
                                         * Referenced by: '<S8>/MATLAB Function'
                                         */
  real_T Controller_SLACK_SPEED;       /* Mask Parameter: Controller_SLACK_SPEED
                                        * Referenced by: '<S1>/Controller'
                                        */
  real_T Controller_SWING_MAX_ANGLE;   /* Mask Parameter: Controller_SWING_MAX_ANGLE
                                        * Referenced by: '<S1>/Controller'
                                        */
  real_T uOrderTD_T;                   /* Mask Parameter: uOrderTD_T
                                        * Referenced by:
                                        *   '<S22>/Gain1'
                                        *   '<S22>/Gain2'
                                        */
  real_T uOrderTD_T1;                  /* Mask Parameter: uOrderTD_T1
                                        * Referenced by:
                                        *   '<S37>/Gain1'
                                        *   '<S37>/Gain2'
                                        *   '<S37>/Gain4'
                                        */
  real_T uOrderTD_T1_l;                /* Mask Parameter: uOrderTD_T1_l
                                        * Referenced by:
                                        *   '<S23>/Gain1'
                                        *   '<S23>/Gain2'
                                        *   '<S23>/Gain4'
                                        */
  real_T uOrderTD_T2;                  /* Mask Parameter: uOrderTD_T2
                                        * Referenced by:
                                        *   '<S37>/Gain1'
                                        *   '<S37>/Gain2'
                                        *   '<S37>/Gain4'
                                        */
  real_T uOrderTD_T2_h;                /* Mask Parameter: uOrderTD_T2_h
                                        * Referenced by:
                                        *   '<S23>/Gain1'
                                        *   '<S23>/Gain2'
                                        *   '<S23>/Gain4'
                                        */
  real_T uOrderTD_Ts;                  /* Mask Parameter: uOrderTD_Ts
                                        * Referenced by:
                                        *   '<S37>/Gain'
                                        *   '<S37>/Gain3'
                                        */
  real_T uOrderTD_Ts_n;                /* Mask Parameter: uOrderTD_Ts_n
                                        * Referenced by:
                                        *   '<S23>/Gain'
                                        *   '<S23>/Gain3'
                                        */
  real_T uOrderTD_Ts_g;                /* Mask Parameter: uOrderTD_Ts_g
                                        * Referenced by: '<S22>/Gain'
                                        */
  real_T Dataprocess_load_vol_gain;    /* Mask Parameter: Dataprocess_load_vol_gain
                                        * Referenced by: '<S21>/Data process'
                                        */
  real_T Dataprocess_load_vol_offset;  /* Mask Parameter: Dataprocess_load_vol_offset
                                        * Referenced by: '<S21>/Data process'
                                        */
  int8_T Controller_run_mode;          /* Mask Parameter: Controller_run_mode
                                        * Referenced by: '<S1>/Controller'
                                        */
  boolean_T StateMachine_BT_CALIB;     /* Mask Parameter: StateMachine_BT_CALIB
                                        * Referenced by: '<S5>/State Machine'
                                        */
  boolean_T StateMachine_BT_ERROR;     /* Mask Parameter: StateMachine_BT_ERROR
                                        * Referenced by: '<S5>/State Machine'
                                        */
  boolean_T StateMachine_BT_IDLE;      /* Mask Parameter: StateMachine_BT_IDLE
                                        * Referenced by: '<S5>/State Machine'
                                        */
  boolean_T Dataprocess_BT_RESET_ANKLE;/* Mask Parameter: Dataprocess_BT_RESET_ANKLE
                                        * Referenced by: '<S19>/Data process'
                                        */
  boolean_T Dataprocess1_BT_RESET_MOTOR;/* Mask Parameter: Dataprocess1_BT_RESET_MOTOR
                                         * Referenced by: '<S19>/Data process1'
                                         */
  boolean_T Dataprocess_BT_RESET_TORQUE;/* Mask Parameter: Dataprocess_BT_RESET_TORQUE
                                         * Referenced by: '<S21>/Data process'
                                         */
  boolean_T StateMachine_BT_RUN;       /* Mask Parameter: StateMachine_BT_RUN
                                        * Referenced by: '<S5>/State Machine'
                                        */
  boolean_T StateMachine_BT_SLACK;     /* Mask Parameter: StateMachine_BT_SLACK
                                        * Referenced by: '<S5>/State Machine'
                                        */
  real_T Gain2_Gain;                   /* Expression: 1/360
                                        * Referenced by: '<S8>/Gain2'
                                        */
  real_T Gain1_Gain;                   /* Expression: 1/20
                                        * Referenced by: '<S8>/Gain1'
                                        */
  real_T Gain_Gain;                    /* Expression: 10
                                        * Referenced by: '<S21>/Gain'
                                        */
  real_T u4low2_NumCoef[4];            /* Expression: [0.219606211225382   0.658818633676145   0.658818633676145   0.219606211225382]*1e-3
                                        * Referenced by: '<S21>/0.4low2'
                                        */
  real_T u4low2_DenCoef[4];            /* Expression: [1.000000000000000  -2.748835809214676   2.528231219142559  -0.777638560238080]
                                        * Referenced by: '<S21>/0.4low2'
                                        */
  real_T u4low2_InitialStates;         /* Expression: 0
                                        * Referenced by: '<S21>/0.4low2'
                                        */
  real_T UnitDelay1_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S37>/Unit Delay1'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S37>/Unit Delay'
                                        */
  real_T UnitDelay2_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S37>/Unit Delay2'
                                        */
  real_T RT1_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<S4>/RT1'
                                        */
  real_T Gain_Gain_c;                  /* Expression: 360/1250
                                        * Referenced by: '<S19>/Gain'
                                        */
  real_T Gain1_Gain_e;                 /* Expression: 360/1250
                                        * Referenced by: '<S19>/Gain1'
                                        */
  real_T RT2_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<S4>/RT2'
                                        */
  real_T Gain2_Gain_e;                 /* Expression: -360.0/1250.0
                                        * Referenced by: '<S19>/Gain2'
                                        */
  real_T Gain3_Gain;                   /* Expression: -360.0/1250.0
                                        * Referenced by: '<S19>/Gain3'
                                        */
  real_T UnitDelay1_InitialCondition_j;/* Expression: 0
                                        * Referenced by: '<S23>/Unit Delay1'
                                        */
  real_T UnitDelay_InitialCondition_j; /* Expression: 0
                                        * Referenced by: '<S23>/Unit Delay'
                                        */
  real_T UnitDelay2_InitialCondition_h;/* Expression: 0
                                        * Referenced by: '<S23>/Unit Delay2'
                                        */
  real_T RT3_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<S4>/RT3'
                                        */
  real_T RT1_InitialCondition_i;       /* Expression: 0
                                        * Referenced by: '<S5>/RT1'
                                        */
  real_T Kp_Value;                     /* Expression: 0
                                        * Referenced by: '<S12>/Kp'
                                        */
  real_T Kd_Value;                     /* Expression: 0
                                        * Referenced by: '<S12>/Kd'
                                        */
  real_T Kl_Value;                     /* Expression: 0
                                        * Referenced by: '<S12>/Kl'
                                        */
  real_T Ks_Value;                     /* Expression: 2
                                        * Referenced by: '<S12>/Ks'
                                        */
  real_T Ko_Value;                     /* Expression: 0
                                        * Referenced by: '<S12>/Ko'
                                        */
  real_T RT1_InitialCondition_e;       /* Expression: 0
                                        * Referenced by: '<S2>/RT1'
                                        */
  real_T peak_torque_Value;            /* Expression: 30
                                        * Referenced by: '<S13>/peak_torque'
                                        */
  real_T rise_time_Value;              /* Expression: 28
                                        * Referenced by: '<S13>/rise_time'
                                        */
  real_T peak_time_Value;              /* Expression: 50
                                        * Referenced by: '<S13>/peak_time'
                                        */
  real_T fall_time_Value;              /* Expression: 15
                                        * Referenced by: '<S13>/fall_time'
                                        */
  real_T RT2_InitialCondition_e;       /* Expression: 0
                                        * Referenced by: '<S2>/RT2'
                                        */
  real_T TSamp_WtEt;                   /* Computed Parameter: TSamp_WtEt
                                        * Referenced by: '<S26>/TSamp'
                                        */
  real_T u4low1_NumCoef[4];            /* Expression: [   0.001567010350588   0.004701031051765   0.004701031051765   0.001567010350588
                                          ]
                                        * Referenced by: '<S19>/0.4low1'
                                        */
  real_T u4low1_DenCoef[4];            /* Expression: [   1.000000000000000  -2.498608344691178   2.115254127003159  -0.604109699507275
                                          ]
                                        * Referenced by: '<S19>/0.4low1'
                                        */
  real_T u4low1_InitialStates;         /* Expression: 0
                                        * Referenced by: '<S19>/0.4low1'
                                        */
  real_T UnitDelay_InitialCondition_g; /* Expression: 0
                                        * Referenced by: '<S22>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition_i;/* Expression: 0
                                        * Referenced by: '<S22>/Unit Delay1'
                                        */
  real_T DataStoreMemory_InitialValue[4400];/* Expression: zeros(4,1100)
                                             * Referenced by: '<Root>/Data Store Memory'
                                             */
  boolean_T VCC1_Value;                /* Computed Parameter: VCC1_Value
                                        * Referenced by: '<S19>/VCC1'
                                        */
  boolean_T VCC3_Value;                /* Computed Parameter: VCC3_Value
                                        * Referenced by: '<S19>/VCC3'
                                        */
  boolean_T Constant_Value;            /* Computed Parameter: Constant_Value
                                        * Referenced by: '<S20>/Constant'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_Human_in_Loop_T {
  const char_T *errorStatus;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T taskTime0;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    boolean_T stopRequestedFlag;
  } Timing;
};

/* Block parameters (auto storage) */
extern P_Human_in_Loop_T Human_in_Loop_P;

/* Block signals (auto storage) */
extern B_Human_in_Loop_T Human_in_Loop_B;

/* Block states (auto storage) */
extern DW_Human_in_Loop_T Human_in_Loop_DW;

/* Model entry point functions */
extern void Human_in_Loop_initialize(void);
extern void Human_in_Loop_output(void);
extern void Human_in_Loop_update(void);
extern void Human_in_Loop_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Human_in_Loop_T *const Human_in_Loop_M;

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
 * '<Root>' : 'Human_in_Loop'
 * '<S1>'   : 'Human_in_Loop/Control Module'
 * '<S2>'   : 'Human_in_Loop/Parameter Module'
 * '<S3>'   : 'Human_in_Loop/RTI Data'
 * '<S4>'   : 'Human_in_Loop/Sensor Data'
 * '<S5>'   : 'Human_in_Loop/State Module'
 * '<S6>'   : 'Human_in_Loop/Timer Interrupt'
 * '<S7>'   : 'Human_in_Loop/Control Module/Controller'
 * '<S8>'   : 'Human_in_Loop/Control Module/Motor'
 * '<S9>'   : 'Human_in_Loop/Control Module/Torque track'
 * '<S10>'  : 'Human_in_Loop/Control Module/Motor/DAC_CLASS1_BL1'
 * '<S11>'  : 'Human_in_Loop/Control Module/Motor/MATLAB Function'
 * '<S12>'  : 'Human_in_Loop/Parameter Module/Control Parameter'
 * '<S13>'  : 'Human_in_Loop/Parameter Module/Torque Parameter'
 * '<S14>'  : 'Human_in_Loop/Parameter Module/Control Parameter/Mux1'
 * '<S15>'  : 'Human_in_Loop/Parameter Module/Torque Parameter/Mux1'
 * '<S16>'  : 'Human_in_Loop/RTI Data/RTI Data Store'
 * '<S17>'  : 'Human_in_Loop/RTI Data/RTI Data Store/RTI Data Store'
 * '<S18>'  : 'Human_in_Loop/RTI Data/RTI Data Store/RTI Data Store/RTI Data Store'
 * '<S19>'  : 'Human_in_Loop/Sensor Data/Encoder module'
 * '<S20>'  : 'Human_in_Loop/Sensor Data/FootSwitch module'
 * '<S21>'  : 'Human_in_Loop/Sensor Data/Torque module'
 * '<S22>'  : 'Human_in_Loop/Sensor Data/Encoder module/1-Order TD'
 * '<S23>'  : 'Human_in_Loop/Sensor Data/Encoder module/2-Order TD'
 * '<S24>'  : 'Human_in_Loop/Sensor Data/Encoder module/Data process'
 * '<S25>'  : 'Human_in_Loop/Sensor Data/Encoder module/Data process1'
 * '<S26>'  : 'Human_in_Loop/Sensor Data/Encoder module/Discrete Derivative'
 * '<S27>'  : 'Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL1'
 * '<S28>'  : 'Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL3'
 * '<S29>'  : 'Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup1'
 * '<S30>'  : 'Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup3'
 * '<S31>'  : 'Human_in_Loop/Sensor Data/Encoder module/Mux'
 * '<S32>'  : 'Human_in_Loop/Sensor Data/Encoder module/Mux1'
 * '<S33>'  : 'Human_in_Loop/Sensor Data/Encoder module/Mux2'
 * '<S34>'  : 'Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_IN_BL1'
 * '<S35>'  : 'Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_OUT_BL1'
 * '<S36>'  : 'Human_in_Loop/Sensor Data/FootSwitch module/FootSwitch Filter'
 * '<S37>'  : 'Human_in_Loop/Sensor Data/Torque module/2-Order TD'
 * '<S38>'  : 'Human_in_Loop/Sensor Data/Torque module/ADC_CLASS1_BL6'
 * '<S39>'  : 'Human_in_Loop/Sensor Data/Torque module/Data process'
 * '<S40>'  : 'Human_in_Loop/Sensor Data/Torque module/MATLAB Function'
 * '<S41>'  : 'Human_in_Loop/Sensor Data/Torque module/Mux'
 * '<S42>'  : 'Human_in_Loop/State Module/Mux1'
 * '<S43>'  : 'Human_in_Loop/State Module/State Machine'
 */
#endif                                 /* RTW_HEADER_Human_in_Loop_h_ */
