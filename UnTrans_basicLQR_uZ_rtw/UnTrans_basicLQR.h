/*
 * UnTrans_basicLQR.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "UnTrans_basicLQR".
 *
 * Model version              : 1.310
 * Simulink Coder version : 9.8 (R2022b) 13-May-2022
 * C source code generated on : Tue Feb  3 15:23:32 2026
 *
 * Target selection: UnTransZynq.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Emulation hardware selection:
 *    Differs from embedded hardware (MATLAB Host)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_UnTrans_basicLQR_h_
#define RTW_HEADER_UnTrans_basicLQR_h_
#ifndef UnTrans_basicLQR_COMMON_INCLUDES_
#define UnTrans_basicLQR_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rt_logging.h"
#include "dt_info.h"
#include "ext_work.h"
#endif                                 /* UnTrans_basicLQR_COMMON_INCLUDES_ */

#include "UnTrans_basicLQR_types.h"
#include <stddef.h>
#include "rt_nonfinite.h"
#include <string.h>
#include "rtGetInf.h"
#include <float.h>

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
#define rtmGetFinalTime(rtm)           ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetRTWExtModeInfo
#define rtmGetRTWExtModeInfo(rtm)      ((rtm)->extModeInfo)
#endif

#ifndef rtmGetRTWLogInfo
#define rtmGetRTWLogInfo(rtm)          ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetSampleHitArray
#define rtmGetSampleHitArray(rtm)      ((rtm)->Timing.sampleHitArray)
#endif

#ifndef rtmGetStepSize
#define rtmGetStepSize(rtm)            ((rtm)->Timing.stepSize)
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGet_TimeOfLastOutput
#define rtmGet_TimeOfLastOutput(rtm)   ((rtm)->Timing.timeOfLastOutput)
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

#ifndef rtmGetTStart
#define rtmGetTStart(rtm)              ((rtm)->Timing.tStart)
#endif

#ifndef rtmGetTimeOfLastOutput
#define rtmGetTimeOfLastOutput(rtm)    ((rtm)->Timing.timeOfLastOutput)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T Encoders[2];                  /* '<S8>/S-Function3' */
  real_T IMU[6];                       /* '<S8>/S-Function3' */
  real_T PWMstatus;                    /* '<S8>/S-Function3' */
  real_T IMUstatus;                    /* '<S8>/S-Function3' */
  real_T theta;                        /* '<S6>/Sum1' */
  real_T DiscreteFilter;               /* '<S6>/Discrete Filter' */
  real_T fi;                           /* '<S6>/Pulse2SI2' */
  real_T DiscreteFilter1;              /* '<S6>/Discrete Filter1' */
  real_T Sum;                          /* '<S3>/Sum' */
  real_T Sum_i;                        /* '<S2>/Sum' */
  real_T ResetEncodersWatchdogandIntegra;
                              /* '<S4>/Reset Encoders Watchdog and Integrals' */
  real_T PWMperiodsec;                 /* '<S8>/Constant1' */
  real_T Gain;                         /* '<S8>/Gain' */
  real_T Gain2;                        /* '<S8>/Gain2' */
  real_T Sumi;                         /* '<S9>/Sumi' */
  real_T SFunction;                    /* '<S15>/S-Function' */
  real_T uDLookupTable;                /* '<S15>/1-D Lookup Table' */
  real_T SFunction2;                   /* '<S15>/S-Function2' */
  real_T Constant2;                    /* '<S15>/Constant2' */
  real_T SFunction3[9];                /* '<S15>/S-Function3' */
  real_T MinMax;                       /* '<S16>/MinMax' */
  real_T MinMax1;                      /* '<S16>/MinMax1' */
  real_T psi;                          /* '<S6>/Kalman' */
  real_T gyroOut;                      /* '<S6>/Kalman' */
} B_UnTrans_basicLQR_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T DiscreteFilter_states;        /* '<S6>/Discrete Filter' */
  real_T DiscreteFilter1_states;       /* '<S6>/Discrete Filter1' */
  real_T Sumi_DSTATE;                  /* '<S9>/Sumi' */
  real_T Sumi1_DSTATE;                 /* '<S15>/Sumi1' */
  real_T UnitDelay_DSTATE;             /* '<S16>/Unit Delay' */
  real_T UnitDelay1_DSTATE;            /* '<S16>/Unit Delay1' */
  real_T Memory2_PreviousInput;        /* '<S12>/Memory2' */
  real_T Memory2_PreviousInput_j;      /* '<S13>/Memory2' */
  real_T Memory1_PreviousInput;        /* '<S15>/Memory1' */
  real_T Memory_PreviousInput;         /* '<S15>/Memory' */
  real_T H[3];                         /* '<S6>/Kalman' */
  real_T Q[9];                         /* '<S6>/Kalman' */
  real_T R;                            /* '<S6>/Kalman' */
  real_T P[9];                         /* '<S6>/Kalman' */
  real_T I[9];                         /* '<S6>/Kalman' */
  real_T x_est[3];                     /* '<S6>/Kalman' */
  real_T q;                            /* '<S6>/Kalman' */
  struct {
    void *LoggedData[2];
  } Ctrl_PWORK;                        /* '<Root>/Ctrl' */

  struct {
    void *LoggedData[6];
  } State_PWORK;                       /* '<Root>/State' */

  int8_T MinMaxPeriod_SubsysRanBC;     /* '<S15>/Min//Max Period' */
  boolean_T Relay_Mode;                /* '<S10>/Relay' */
  boolean_T P_not_empty;               /* '<S6>/Kalman' */
} DW_UnTrans_basicLQR_T;

/* Parameters (default storage) */
struct P_UnTrans_basicLQR_T_ {
  real_T SFunction3_P1_Size[2];        /* Computed Parameter: SFunction3_P1_Size
                                        * Referenced by:
                                        *   '<S8>/S-Function3'
                                        *   '<S15>/S-Function'
                                        *   '<S15>/S-Function1'
                                        *   '<S15>/S-Function2'
                                        *   '<S15>/S-Function3'
                                        */
  real_T Ts;                           /* Variable: Ts
                                        * Referenced by:
                                        *   '<S6>/Constant'
                                        *   '<S8>/S-Function3'
                                        *   '<S12>/Constant'
                                        *   '<S13>/Constant'
                                        *   '<S15>/S-Function'
                                        *   '<S15>/S-Function1'
                                        *   '<S15>/S-Function2'
                                        *   '<S15>/S-Function3'
                                        */
  real_T Translation_gain;             /* Mask Parameter: Translation_gain
                                        * Referenced by: '<S7>/Slider Gain'
                                        */
  real_T Rotation_gain;                /* Mask Parameter: Rotation_gain
                                        * Referenced by: '<S5>/Slider Gain'
                                        */
  real_T CoulombViscousFriction2_gain;
                                 /* Mask Parameter: CoulombViscousFriction2_gain
                                  * Referenced by: '<S3>/Gain'
                                  */
  real_T CoulombViscousFriction1_gain;
                                 /* Mask Parameter: CoulombViscousFriction1_gain
                                  * Referenced by: '<S2>/Gain'
                                  */
  real_T CoulombViscousFriction2_offset;
                               /* Mask Parameter: CoulombViscousFriction2_offset
                                * Referenced by: '<S3>/Gain1'
                                */
  real_T CoulombViscousFriction1_offset;
                               /* Mask Parameter: CoulombViscousFriction1_offset
                                * Referenced by: '<S2>/Gain1'
                                */
  real_T NoCtrl_Value[2];              /* Expression: [0; 0]
                                        * Referenced by: '<S1>/NoCtrl'
                                        */
  real_T deg2rad_Gain;                 /* Expression: pi/180
                                        * Referenced by: '<S1>/deg2rad'
                                        */
  real_T Uprightanglecorrection_Value; /* Expression: 0*pi/180
                                        * Referenced by: '<Root>/Upright angle correction'
                                        */
  real_T u_Value;                      /* Expression: 0
                                        * Referenced by: '<S1>/0'
                                        */
  real_T m2rad_Gain;                   /* Expression: 13.3333
                                        * Referenced by: '<S1>/m2rad'
                                        */
  real_T LQR_Gain[12];                 /* Expression: [K  ;K(1:4) -K(5:6)]
                                        * Referenced by: '<S11>/LQR'
                                        */
  real_T Normal_Value;                 /* Expression: 0
                                        * Referenced by: '<S4>/Normal'
                                        */
  real_T Switch_Threshold;             /* Expression: 0.01
                                        * Referenced by: '<S4>/Switch'
                                        */
  real_T Mu2_Gain;                     /* Expression: 1000
                                        * Referenced by: '<S16>/Mu2'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 9999
                                        * Referenced by: '<S16>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition;  /* Expression: -9999
                                        * Referenced by: '<S16>/Unit Delay1'
                                        */
  real_T raw2g_Gain;                   /* Expression: 4/32767
                                        * Referenced by: '<S8>/raw2g'
                                        */
  real_T Gain4_Gain;                   /* Expression: -1
                                        * Referenced by: '<S8>/Gain4'
                                        */
  real_T raw2dps_Gain;                 /* Expression: 250/32767
                                        * Referenced by: '<S8>/raw2dps'
                                        */
  real_T UnitConversion_Gain;          /* Expression: pi/180
                                        * Referenced by: '<S6>/UnitConversion'
                                        */
  real_T Relay_OnVal;                  /* Expression: 12*pi/180
                                        * Referenced by: '<S10>/Relay'
                                        */
  real_T Relay_OffVal;                 /* Expression: 2*pi/180
                                        * Referenced by: '<S10>/Relay'
                                        */
  real_T Relay_YOn;                    /* Expression: 1
                                        * Referenced by: '<S10>/Relay'
                                        */
  real_T Relay_YOff;                   /* Expression: 0
                                        * Referenced by: '<S10>/Relay'
                                        */
  real_T Gain3_Gain;                   /* Expression: 1
                                        * Referenced by: '<S8>/Gain3'
                                        */
  real_T Gain1_Gain;                   /* Expression: -1
                                        * Referenced by: '<S8>/Gain1'
                                        */
  real_T Pulse2SI1_Gain;               /* Expression: pi/16000
                                        * Referenced by: '<S6>/Pulse2SI1'
                                        */
  real_T Memory2_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S12>/Memory2'
                                        */
  real_T Saturation_UpperSat;          /* Expression: 1
                                        * Referenced by: '<S12>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: 0.001
                                        * Referenced by: '<S12>/Saturation'
                                        */
  real_T DiscreteFilter_NumCoef[2];    /* Expression: [1 1]
                                        * Referenced by: '<S6>/Discrete Filter'
                                        */
  real_T DiscreteFilter_DenCoef[2];    /* Expression: [5 -3]
                                        * Referenced by: '<S6>/Discrete Filter'
                                        */
  real_T DiscreteFilter_InitialStates; /* Expression: 0
                                        * Referenced by: '<S6>/Discrete Filter'
                                        */
  real_T Pulse2SI2_Gain;               /* Expression: pi/16000*15/40
                                        * Referenced by: '<S6>/Pulse2SI2'
                                        */
  real_T Memory2_InitialCondition_i;   /* Expression: 0
                                        * Referenced by: '<S13>/Memory2'
                                        */
  real_T Saturation_UpperSat_i;        /* Expression: 1
                                        * Referenced by: '<S13>/Saturation'
                                        */
  real_T Saturation_LowerSat_f;        /* Expression: 0.001
                                        * Referenced by: '<S13>/Saturation'
                                        */
  real_T DiscreteFilter1_NumCoef[2];   /* Expression: [1 1]
                                        * Referenced by: '<S6>/Discrete Filter1'
                                        */
  real_T DiscreteFilter1_DenCoef[2];   /* Expression: [5 -3]
                                        * Referenced by: '<S6>/Discrete Filter1'
                                        */
  real_T DiscreteFilter1_InitialStates;/* Expression: 0
                                        * Referenced by: '<S6>/Discrete Filter1'
                                        */
  real_T Stepm_Value;                  /* Expression: 0.5
                                        * Referenced by: '<Root>/Step [m]'
                                        */
  real_T Stepdeg_Value;                /* Expression: 90
                                        * Referenced by: '<Root>/Step [deg]'
                                        */
  real_T Switch_Threshold_h;           /* Expression: 1
                                        * Referenced by: '<S1>/Switch'
                                        */
  real_T Gain_Gain;                    /* Expression: 1
                                        * Referenced by: '<Root>/Gain'
                                        */
  real_T Saturation1_UpperSat;         /* Expression: 0.6
                                        * Referenced by: '<Root>/Saturation1'
                                        */
  real_T Saturation1_LowerSat;         /* Expression: -0.6
                                        * Referenced by: '<Root>/Saturation1'
                                        */
  real_T Saturation3_UpperSat;         /* Expression: 0.6
                                        * Referenced by: '<Root>/Saturation3'
                                        */
  real_T Saturation3_LowerSat;         /* Expression: -0.6
                                        * Referenced by: '<Root>/Saturation3'
                                        */
  real_T Reset_Value;                  /* Expression: 1
                                        * Referenced by: '<S4>/Reset'
                                        */
  real_T Constant1_Value;              /* Expression: 0.0001
                                        * Referenced by: '<S8>/Constant1'
                                        */
  real_T Gain_Gain_n;                  /* Expression: 1
                                        * Referenced by: '<S8>/Gain'
                                        */
  real_T Gain2_Gain;                   /* Expression: -1
                                        * Referenced by: '<S8>/Gain2'
                                        */
  real_T Sumi_gainval;                 /* Computed Parameter: Sumi_gainval
                                        * Referenced by: '<S9>/Sumi'
                                        */
  real_T Sumi_IC;                      /* Expression: 0
                                        * Referenced by: '<S9>/Sumi'
                                        */
  real_T Constant4_Value;              /* Expression: 1
                                        * Referenced by: '<S9>/Constant4'
                                        */
  real_T Sumi1_gainval;                /* Computed Parameter: Sumi1_gainval
                                        * Referenced by: '<S15>/Sumi1'
                                        */
  real_T Sumi1_IC;                     /* Expression: 0
                                        * Referenced by: '<S15>/Sumi1'
                                        */
  real_T uDLookupTable_tableData[8];   /* Expression: [1,1,0,0,1,1,0,0]
                                        * Referenced by: '<S15>/1-D Lookup Table'
                                        */
  real_T uDLookupTable_bp01Data[8];
                               /* Expression: [0,0.05,0.06,0.2,0.21,0.25,0.26,1]
                                * Referenced by: '<S15>/1-D Lookup Table'
                                */
  real_T Memory1_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S15>/Memory1'
                                        */
  real_T Constant1_Value_f;            /* Expression: 1
                                        * Referenced by: '<S15>/Constant1'
                                        */
  real_T Memory_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S15>/Memory'
                                        */
  real_T Mu1_Gain;                     /* Expression: 1/100e6
                                        * Referenced by: '<S15>/Mu1'
                                        */
  real_T Mu2_Gain_p;                   /* Expression: 1000
                                        * Referenced by: '<S15>/Mu2'
                                        */
  real_T Constant2_Value;              /* Expression: 0
                                        * Referenced by: '<S15>/Constant2'
                                        */
  uint8_T ResetEncodersWatchdogandIntegra;
                          /* Computed Parameter: ResetEncodersWatchdogandIntegra
                           * Referenced by: '<S4>/Reset Encoders Watchdog and Integrals'
                           */
};

/* Real-time Model Data Structure */
struct tag_RTM_UnTrans_basicLQR_T {
  struct SimStruct_tag * *childSfunctions;
  const char_T *errorStatus;
  SS_SimMode simMode;
  RTWLogInfo *rtwLogInfo;
  RTWExtModeInfo *extModeInfo;
  RTWSolverInfo solverInfo;
  RTWSolverInfo *solverInfoPtr;
  void *sfcnInfo;

  /*
   * NonInlinedSFcns:
   * The following substructure contains information regarding
   * non-inlined s-functions used in the model.
   */
  struct {
    RTWSfcnInfo sfcnInfo;
    time_T *taskTimePtrs[1];
    SimStruct childSFunctions[5];
    SimStruct *childSFunctionPtrs[5];
    struct _ssBlkInfo2 blkInfo2[5];
    struct _ssSFcnModelMethods2 methods2[5];
    struct _ssSFcnModelMethods3 methods3[5];
    struct _ssSFcnModelMethods4 methods4[5];
    struct _ssStatesInfo2 statesInfo2[5];
    ssPeriodicStatesInfo periodicStatesInfo[5];
    struct _ssPortInfo2 inputOutputPortInfo2[5];
    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[3];
      struct _ssInPortUnit inputPortUnits[3];
      struct _ssInPortCoSimAttribute inputPortCoSimAttribute[3];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      real_T const *UPtrs2[2];
      struct _ssPortOutputs outputPortInfo[4];
      struct _ssOutPortUnit outputPortUnits[4];
      struct _ssOutPortCoSimAttribute outputPortCoSimAttribute[4];
      uint_T attribs[1];
      mxArray *params[1];
    } Sfcn0;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[1];
      struct _ssOutPortUnit outputPortUnits[1];
      struct _ssOutPortCoSimAttribute outputPortCoSimAttribute[1];
      uint_T attribs[1];
      mxArray *params[1];
    } Sfcn1;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      struct _ssInPortUnit inputPortUnits[1];
      struct _ssInPortCoSimAttribute inputPortCoSimAttribute[1];
      real_T const *UPtrs0[1];
      uint_T attribs[1];
      mxArray *params[1];
    } Sfcn2;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortOutputs outputPortInfo[1];
      struct _ssOutPortUnit outputPortUnits[1];
      struct _ssOutPortCoSimAttribute outputPortCoSimAttribute[1];
      uint_T attribs[1];
      mxArray *params[1];
    } Sfcn3;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      struct _ssInPortUnit inputPortUnits[1];
      struct _ssInPortCoSimAttribute inputPortCoSimAttribute[1];
      real_T const *UPtrs0[1];
      struct _ssPortOutputs outputPortInfo[1];
      struct _ssOutPortUnit outputPortUnits[1];
      struct _ssOutPortCoSimAttribute outputPortCoSimAttribute[1];
      uint_T attribs[1];
      mxArray *params[1];
    } Sfcn4;
  } NonInlinedSFcns;

  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T checksums[4];
    uint32_T options;
    int_T numContStates;
    int_T numU;
    int_T numY;
    int_T numSampTimes;
    int_T numBlocks;
    int_T numBlockIO;
    int_T numBlockPrms;
    int_T numDwork;
    int_T numSFcnPrms;
    int_T numSFcns;
    int_T numIports;
    int_T numOports;
    int_T numNonSampZCs;
    int_T sysDirFeedThru;
    int_T rtwGenSfcn;
  } Sizes;

  /*
   * SpecialInfo:
   * The following substructure contains special information
   * related to other components that are dependent on RTW.
   */
  struct {
    const void *mappingInfo;
  } SpecialInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T stepSize;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    time_T tStart;
    time_T tFinal;
    time_T timeOfLastOutput;
    boolean_T stopRequestedFlag;
    time_T *sampleTimes;
    time_T *offsetTimes;
    int_T *sampleTimeTaskIDPtr;
    int_T *sampleHits;
    int_T *perTaskSampleHits;
    time_T *t;
    time_T sampleTimesArray[1];
    time_T offsetTimesArray[1];
    int_T sampleTimeTaskIDArray[1];
    int_T sampleHitArray[1];
    int_T perTaskSampleHitsArray[1];
    time_T tArray[1];
  } Timing;
};

/* Block parameters (default storage) */
extern P_UnTrans_basicLQR_T UnTrans_basicLQR_P;

/* Block signals (default storage) */
extern B_UnTrans_basicLQR_T UnTrans_basicLQR_B;

/* Block states (default storage) */
extern DW_UnTrans_basicLQR_T UnTrans_basicLQR_DW;

/* Model entry point functions */
extern void UnTrans_basicLQR_initialize(void);
extern void UnTrans_basicLQR_step(void);
extern void UnTrans_basicLQR_terminate(void);

/* Real-time Model object */
extern RT_MODEL_UnTrans_basicLQR_T *const UnTrans_basicLQR_M;

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
 * '<Root>' : 'UnTrans_basicLQR'
 * '<S1>'   : 'UnTrans_basicLQR/Controller'
 * '<S2>'   : 'UnTrans_basicLQR/Coulomb & Viscous Friction1'
 * '<S3>'   : 'UnTrans_basicLQR/Coulomb & Viscous Friction2'
 * '<S4>'   : 'UnTrans_basicLQR/Reset'
 * '<S5>'   : 'UnTrans_basicLQR/Rotation'
 * '<S6>'   : 'UnTrans_basicLQR/StateObserver'
 * '<S7>'   : 'UnTrans_basicLQR/Translation'
 * '<S8>'   : 'UnTrans_basicLQR/UnTrans'
 * '<S9>'   : 'UnTrans_basicLQR/microZed '
 * '<S10>'  : 'UnTrans_basicLQR/Controller/Safety'
 * '<S11>'  : 'UnTrans_basicLQR/Controller/Subsystem'
 * '<S12>'  : 'UnTrans_basicLQR/StateObserver/Derivative'
 * '<S13>'  : 'UnTrans_basicLQR/StateObserver/Derivative1'
 * '<S14>'  : 'UnTrans_basicLQR/StateObserver/Kalman'
 * '<S15>'  : 'UnTrans_basicLQR/microZed /MicroZed'
 * '<S16>'  : 'UnTrans_basicLQR/microZed /MicroZed/Min//Max Period'
 */
#endif                                 /* RTW_HEADER_UnTrans_basicLQR_h_ */
