/*
 * UnTrans_basicLQR_data.c
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

#include "UnTrans_basicLQR.h"

/* Block parameters (default storage) */
P_UnTrans_basicLQR_T UnTrans_basicLQR_P = {
  /* Computed Parameter: SFunction3_P1_Size
   * Referenced by:
   *   '<S8>/S-Function3'
   *   '<S15>/S-Function'
   *   '<S15>/S-Function1'
   *   '<S15>/S-Function2'
   *   '<S15>/S-Function3'
   */
  { 1.0, 1.0 },

  /* Variable: Ts
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
  0.01,

  /* Mask Parameter: Translation_gain
   * Referenced by: '<S7>/Slider Gain'
   */
  0.0,

  /* Mask Parameter: Rotation_gain
   * Referenced by: '<S5>/Slider Gain'
   */
  0.0,

  /* Mask Parameter: CoulombViscousFriction2_gain
   * Referenced by: '<S3>/Gain'
   */
  1.0,

  /* Mask Parameter: CoulombViscousFriction1_gain
   * Referenced by: '<S2>/Gain'
   */
  1.0,

  /* Mask Parameter: CoulombViscousFriction2_offset
   * Referenced by: '<S3>/Gain1'
   */
  0.055,

  /* Mask Parameter: CoulombViscousFriction1_offset
   * Referenced by: '<S2>/Gain1'
   */
  0.055,

  /* Expression: [0; 0]
   * Referenced by: '<S1>/NoCtrl'
   */
  { 0.0, 0.0 },

  /* Expression: pi/180
   * Referenced by: '<S1>/deg2rad'
   */
  0.017453292519943295,

  /* Expression: 0*pi/180
   * Referenced by: '<Root>/Upright angle correction'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S1>/0'
   */
  0.0,

  /* Expression: 13.3333
   * Referenced by: '<S1>/m2rad'
   */
  13.3333,

  /* Expression: [K  ;K(1:4) -K(5:6)]
   * Referenced by: '<S11>/LQR'
   */
  { -0.0158, -0.0158, -0.0514, -0.0514, -4.1139, -4.1139, -0.216, -0.216, -0.1,
    0.1, -0.11, 0.11 },

  /* Expression: 0
   * Referenced by: '<S4>/Normal'
   */
  0.0,

  /* Expression: 0.01
   * Referenced by: '<S4>/Switch'
   */
  0.01,

  /* Expression: 1000
   * Referenced by: '<S16>/Mu2'
   */
  1000.0,

  /* Expression: 9999
   * Referenced by: '<S16>/Unit Delay'
   */
  9999.0,

  /* Expression: -9999
   * Referenced by: '<S16>/Unit Delay1'
   */
  -9999.0,

  /* Expression: 4/32767
   * Referenced by: '<S8>/raw2g'
   */
  0.00012207403790398877,

  /* Expression: -1
   * Referenced by: '<S8>/Gain4'
   */
  -1.0,

  /* Expression: 250/32767
   * Referenced by: '<S8>/raw2dps'
   */
  0.0076296273689992981,

  /* Expression: pi/180
   * Referenced by: '<S6>/UnitConversion'
   */
  0.017453292519943295,

  /* Expression: 12*pi/180
   * Referenced by: '<S10>/Relay'
   */
  0.20943951023931953,

  /* Expression: 2*pi/180
   * Referenced by: '<S10>/Relay'
   */
  0.034906585039886591,

  /* Expression: 1
   * Referenced by: '<S10>/Relay'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S10>/Relay'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S8>/Gain3'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S8>/Gain1'
   */
  -1.0,

  /* Expression: pi/16000
   * Referenced by: '<S6>/Pulse2SI1'
   */
  0.00019634954084936208,

  /* Expression: 0
   * Referenced by: '<S12>/Memory2'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S12>/Saturation'
   */
  1.0,

  /* Expression: 0.001
   * Referenced by: '<S12>/Saturation'
   */
  0.001,

  /* Expression: [1 1]
   * Referenced by: '<S6>/Discrete Filter'
   */
  { 1.0, 1.0 },

  /* Expression: [5 -3]
   * Referenced by: '<S6>/Discrete Filter'
   */
  { 5.0, -3.0 },

  /* Expression: 0
   * Referenced by: '<S6>/Discrete Filter'
   */
  0.0,

  /* Expression: pi/16000*15/40
   * Referenced by: '<S6>/Pulse2SI2'
   */
  7.3631077818510777E-5,

  /* Expression: 0
   * Referenced by: '<S13>/Memory2'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S13>/Saturation'
   */
  1.0,

  /* Expression: 0.001
   * Referenced by: '<S13>/Saturation'
   */
  0.001,

  /* Expression: [1 1]
   * Referenced by: '<S6>/Discrete Filter1'
   */
  { 1.0, 1.0 },

  /* Expression: [5 -3]
   * Referenced by: '<S6>/Discrete Filter1'
   */
  { 5.0, -3.0 },

  /* Expression: 0
   * Referenced by: '<S6>/Discrete Filter1'
   */
  0.0,

  /* Expression: 0.5
   * Referenced by: '<Root>/Step [m]'
   */
  0.5,

  /* Expression: 90
   * Referenced by: '<Root>/Step [deg]'
   */
  90.0,

  /* Expression: 1
   * Referenced by: '<S1>/Switch'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<Root>/Gain'
   */
  1.0,

  /* Expression: 0.6
   * Referenced by: '<Root>/Saturation1'
   */
  0.6,

  /* Expression: -0.6
   * Referenced by: '<Root>/Saturation1'
   */
  -0.6,

  /* Expression: 0.6
   * Referenced by: '<Root>/Saturation3'
   */
  0.6,

  /* Expression: -0.6
   * Referenced by: '<Root>/Saturation3'
   */
  -0.6,

  /* Expression: 1
   * Referenced by: '<S4>/Reset'
   */
  1.0,

  /* Expression: 0.0001
   * Referenced by: '<S8>/Constant1'
   */
  0.0001,

  /* Expression: 1
   * Referenced by: '<S8>/Gain'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S8>/Gain2'
   */
  -1.0,

  /* Computed Parameter: Sumi_gainval
   * Referenced by: '<S9>/Sumi'
   */
  0.01,

  /* Expression: 0
   * Referenced by: '<S9>/Sumi'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S9>/Constant4'
   */
  1.0,

  /* Computed Parameter: Sumi1_gainval
   * Referenced by: '<S15>/Sumi1'
   */
  0.01,

  /* Expression: 0
   * Referenced by: '<S15>/Sumi1'
   */
  0.0,

  /* Expression: [1,1,0,0,1,1,0,0]
   * Referenced by: '<S15>/1-D Lookup Table'
   */
  { 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0 },

  /* Expression: [0,0.05,0.06,0.2,0.21,0.25,0.26,1]
   * Referenced by: '<S15>/1-D Lookup Table'
   */
  { 0.0, 0.05, 0.06, 0.2, 0.21, 0.25, 0.26, 1.0 },

  /* Expression: 0
   * Referenced by: '<S15>/Memory1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S15>/Constant1'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S15>/Memory'
   */
  0.0,

  /* Expression: 1/100e6
   * Referenced by: '<S15>/Mu1'
   */
  1.0E-8,

  /* Expression: 1000
   * Referenced by: '<S15>/Mu2'
   */
  1000.0,

  /* Expression: 0
   * Referenced by: '<S15>/Constant2'
   */
  0.0,

  /* Computed Parameter: ResetEncodersWatchdogandIntegra
   * Referenced by: '<S4>/Reset Encoders Watchdog and Integrals'
   */
  0U
};
