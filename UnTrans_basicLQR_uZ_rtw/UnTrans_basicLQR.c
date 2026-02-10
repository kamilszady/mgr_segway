/*
 * UnTrans_basicLQR.c
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
#include "rtwtypes.h"
#include "UnTrans_basicLQR_private.h"
#include <string.h>
#include <math.h>
#include "rt_nonfinite.h"
#include "rt_defines.h"
#include <float.h>
#include "UnTrans_basicLQR_dt.h"

/* Block signals (default storage) */
B_UnTrans_basicLQR_T UnTrans_basicLQR_B;

/* Block states (default storage) */
DW_UnTrans_basicLQR_T UnTrans_basicLQR_DW;

/* Real-time model */
static RT_MODEL_UnTrans_basicLQR_T UnTrans_basicLQR_M_;
RT_MODEL_UnTrans_basicLQR_T *const UnTrans_basicLQR_M = &UnTrans_basicLQR_M_;
real_T look1_binlxpw(real_T u0, const real_T bp0[], const real_T table[],
                     uint32_T maxIndex)
{
  real_T frac;
  real_T yL_0d0;
  uint32_T bpIdx;
  uint32_T iLeft;
  uint32_T iRght;

  /* Column-major Lookup 1-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex]) {
    /* Binary Search */
    bpIdx = maxIndex >> 1U;
    iLeft = 0U;
    iRght = maxIndex;
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex - 1U;
    frac = (u0 - bp0[maxIndex - 1U]) / (bp0[maxIndex] - bp0[maxIndex - 1U]);
  }

  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  yL_0d0 = table[iLeft];
  return (table[iLeft + 1U] - yL_0d0) * frac + yL_0d0;
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  int32_T tmp;
  int32_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      tmp = 1;
    } else {
      tmp = -1;
    }

    if (u1 > 0.0) {
      tmp_0 = 1;
    } else {
      tmp_0 = -1;
    }

    y = atan2(tmp, tmp_0);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

real_T rt_remd_snf(real_T u0, real_T u1)
{
  real_T q;
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = (rtNaN);
  } else if (rtIsInf(u1)) {
    y = u0;
  } else {
    if (u1 < 0.0) {
      q = ceil(u1);
    } else {
      q = floor(u1);
    }

    if ((u1 != 0.0) && (u1 != q)) {
      q = fabs(u0 / u1);
      if (!(fabs(q - floor(q + 0.5)) > DBL_EPSILON * q)) {
        y = 0.0 * u0;
      } else {
        y = fmod(u0, u1);
      }
    } else {
      y = fmod(u0, u1);
    }
  }

  return y;
}

/* Model step function */
void UnTrans_basicLQR_step(void)
{
  /* local block i/o variables */
  real_T rtb_Add2;
  real_T rtb_Mu2;
  real_T A[9];
  real_T A_0[9];
  real_T tmp[6];
  real_T tmp_0[6];
  real_T tmp_1[6];
  real_T K[3];
  real_T dt_0[3];
  real_T rtb_Gain_a[2];
  real_T dt;
  real_T rtb_Fcn1;
  real_T rtb_Sumi1;
  real_T rtb_UnitConversion;
  real_T rtb_UnitDelay1;
  int32_T A_tmp;
  int32_T i;
  int32_T i_0;

  /* Reset subsysRan breadcrumbs */
  srClearBC(UnTrans_basicLQR_DW.MinMaxPeriod_SubsysRanBC);

  /* S-Function (microZedUnTrans): '<S8>/S-Function3' */

  /* Level2 S-Function Block: '<S8>/S-Function3' (microZedUnTrans) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[0];
    sfcnOutputs(rts,0);
  }

  /* Fcn: '<S6>/Fcn1' incorporates:
   *  Gain: '<S8>/Gain4'
   *  Gain: '<S8>/raw2g'
   */
  rtb_Fcn1 = rt_atan2d_snf(UnTrans_basicLQR_P.raw2g_Gain *
    UnTrans_basicLQR_B.IMU[0] * UnTrans_basicLQR_P.Gain4_Gain,
    UnTrans_basicLQR_P.raw2g_Gain * UnTrans_basicLQR_B.IMU[2]);

  /* Gain: '<S6>/UnitConversion' incorporates:
   *  Gain: '<S8>/raw2dps'
   */
  rtb_UnitConversion = UnTrans_basicLQR_P.raw2dps_Gain * UnTrans_basicLQR_B.IMU
    [4] * UnTrans_basicLQR_P.UnitConversion_Gain;

  /* MATLAB Function: '<S6>/Kalman' incorporates:
   *  Constant: '<S6>/Constant'
   */
  /* MATLAB Function 'StateObserver/Kalman': '<S14>:1' */
  if (!UnTrans_basicLQR_DW.P_not_empty) {
    /* '<S14>:1:5' */
    /* '<S14>:1:8' */
    UnTrans_basicLQR_DW.Q[0] = UnTrans_basicLQR_DW.q;
    UnTrans_basicLQR_DW.Q[3] = 0.0;
    UnTrans_basicLQR_DW.Q[6] = 0.0;
    UnTrans_basicLQR_DW.Q[1] = 0.0;
    UnTrans_basicLQR_DW.Q[4] = UnTrans_basicLQR_DW.q;
    UnTrans_basicLQR_DW.Q[7] = 0.0;
    UnTrans_basicLQR_DW.Q[2] = 0.0;
    UnTrans_basicLQR_DW.Q[5] = 0.0;
    UnTrans_basicLQR_DW.Q[8] = UnTrans_basicLQR_DW.q;

    /* '<S14>:1:12' */
    UnTrans_basicLQR_DW.x_est[0] = rtb_Fcn1;
    UnTrans_basicLQR_DW.x_est[1] = 0.0;
    UnTrans_basicLQR_DW.x_est[2] = rtb_UnitConversion;

    /* '<S14>:1:14' */
    UnTrans_basicLQR_DW.P_not_empty = true;

    /* '<S14>:1:15' */
    for (i = 0; i < 9; i++) {
      UnTrans_basicLQR_DW.P[i] = UnTrans_basicLQR_DW.Q[i];
      UnTrans_basicLQR_DW.I[i] = 0.0;
    }

    UnTrans_basicLQR_DW.I[0] = 1.0;
    UnTrans_basicLQR_DW.I[4] = 1.0;
    UnTrans_basicLQR_DW.I[8] = 1.0;
  }

  /* '<S14>:1:18' */
  dt = UnTrans_basicLQR_P.Ts;
  if (UnTrans_basicLQR_P.Ts < 0.002) {
    /* '<S14>:1:19' */
    /* '<S14>:1:20' */
    dt = 0.002;
  } else if (UnTrans_basicLQR_P.Ts > 0.01) {
    /* '<S14>:1:21' */
    /* '<S14>:1:22' */
    dt = 0.01;
  }

  /* '<S14>:1:25' */
  A[0] = 1.0;
  A[3] = 0.0;
  A[6] = -dt;
  A[1] = 0.0;
  A[2] = 0.0;
  A[4] = 0.0;
  A[5] = 0.0;
  A[7] = -1.0;
  A[8] = 1.0;

  /* '<S14>:1:26' */
  /* '<S14>:1:29' */
  for (i = 0; i < 3; i++) {
    K[i] = A[i + 6] * UnTrans_basicLQR_DW.x_est[2] + (0.0 *
      UnTrans_basicLQR_DW.x_est[1] + A[i] * UnTrans_basicLQR_DW.x_est[0]);
  }

  dt_0[0] = dt * rtb_UnitConversion;
  dt_0[1] = rtb_UnitConversion;
  dt_0[2] = 0.0 * rtb_UnitConversion;

  /* '<S14>:1:30' */
  for (i = 0; i < 3; i++) {
    UnTrans_basicLQR_DW.x_est[i] = K[i] + dt_0[i];
    for (i_0 = 0; i_0 < 3; i_0++) {
      A_tmp = 3 * i_0 + i;
      A_0[A_tmp] = 0.0;
      A_0[A_tmp] += UnTrans_basicLQR_DW.P[3 * i_0] * A[i];
      A_0[A_tmp] += UnTrans_basicLQR_DW.P[3 * i_0 + 1] * 0.0;
      A_0[A_tmp] += UnTrans_basicLQR_DW.P[3 * i_0 + 2] * A[i + 6];
    }
  }

  for (i = 0; i < 3; i++) {
    for (i_0 = 0; i_0 < 3; i_0++) {
      A_tmp = 3 * i_0 + i;
      UnTrans_basicLQR_DW.P[A_tmp] = ((A_0[i + 3] * 0.0 + A_0[i] * A[i_0]) +
        A_0[i + 6] * A[i_0 + 6]) + UnTrans_basicLQR_DW.Q[A_tmp];
    }
  }

  /* '<S14>:1:33' */
  rtb_UnitConversion = 0.0;
  for (i = 0; i < 3; i++) {
    rtb_UnitConversion += ((UnTrans_basicLQR_DW.P[3 * i + 1] *
      UnTrans_basicLQR_DW.H[1] + UnTrans_basicLQR_DW.P[3 * i] *
      UnTrans_basicLQR_DW.H[0]) + UnTrans_basicLQR_DW.P[3 * i + 2] *
      UnTrans_basicLQR_DW.H[2]) * UnTrans_basicLQR_DW.H[i];
  }

  rtb_UnitConversion += UnTrans_basicLQR_DW.R;

  /* '<S14>:1:34' */
  rtb_Sumi1 = 0.0;
  for (i = 0; i < 3; i++) {
    K[i] = ((UnTrans_basicLQR_DW.P[i + 3] * UnTrans_basicLQR_DW.H[1] +
             UnTrans_basicLQR_DW.P[i] * UnTrans_basicLQR_DW.H[0]) +
            UnTrans_basicLQR_DW.P[i + 6] * UnTrans_basicLQR_DW.H[2]) /
      rtb_UnitConversion;
    rtb_Sumi1 += UnTrans_basicLQR_DW.H[i] * UnTrans_basicLQR_DW.x_est[i];
  }

  rtb_Fcn1 -= rtb_Sumi1;

  /* '<S14>:1:35' */
  for (i = 0; i < 3; i++) {
    UnTrans_basicLQR_DW.x_est[i] += K[i] * rtb_Fcn1;
    A[3 * i] = UnTrans_basicLQR_DW.I[3 * i] - K[0] * UnTrans_basicLQR_DW.H[i];
    i_0 = 3 * i + 1;
    A[i_0] = UnTrans_basicLQR_DW.I[i_0] - K[1] * UnTrans_basicLQR_DW.H[i];
    i_0 = 3 * i + 2;
    A[i_0] = UnTrans_basicLQR_DW.I[i_0] - K[2] * UnTrans_basicLQR_DW.H[i];
  }

  for (i = 0; i < 3; i++) {
    for (i_0 = 0; i_0 < 3; i_0++) {
      A_tmp = 3 * i + i_0;
      A_0[A_tmp] = 0.0;
      A_0[A_tmp] += UnTrans_basicLQR_DW.P[3 * i] * A[i_0];
      A_0[A_tmp] += UnTrans_basicLQR_DW.P[3 * i + 1] * A[i_0 + 3];
      A_0[A_tmp] += UnTrans_basicLQR_DW.P[3 * i + 2] * A[i_0 + 6];
    }
  }

  memcpy(&UnTrans_basicLQR_DW.P[0], &A_0[0], 9U * sizeof(real_T));

  /* '<S14>:1:37' */
  UnTrans_basicLQR_B.psi = (UnTrans_basicLQR_DW.H[0] *
    UnTrans_basicLQR_DW.x_est[0] + UnTrans_basicLQR_DW.H[1] *
    UnTrans_basicLQR_DW.x_est[1]) + UnTrans_basicLQR_DW.H[2] *
    UnTrans_basicLQR_DW.x_est[2];

  /* '<S14>:1:38' */
  UnTrans_basicLQR_B.gyroOut = UnTrans_basicLQR_DW.x_est[1];

  /* End of MATLAB Function: '<S6>/Kalman' */

  /* Abs: '<S10>/Abs' */
  rtb_UnitConversion = fabs(UnTrans_basicLQR_B.psi);

  /* Relay: '<S10>/Relay' */
  UnTrans_basicLQR_DW.Relay_Mode = ((rtb_UnitConversion >=
    UnTrans_basicLQR_P.Relay_OnVal) || ((!(rtb_UnitConversion <=
    UnTrans_basicLQR_P.Relay_OffVal)) && UnTrans_basicLQR_DW.Relay_Mode));

  /* Gain: '<S8>/Gain3' */
  dt = UnTrans_basicLQR_P.Gain3_Gain * UnTrans_basicLQR_B.Encoders[0];

  /* Gain: '<S8>/Gain1' */
  rtb_UnitDelay1 = UnTrans_basicLQR_P.Gain1_Gain * UnTrans_basicLQR_B.Encoders[1];

  /* Sum: '<S6>/Sum1' incorporates:
   *  Gain: '<S6>/Pulse2SI1'
   *  Sum: '<S6>/Add'
   */
  UnTrans_basicLQR_B.theta = (dt + rtb_UnitDelay1) *
    UnTrans_basicLQR_P.Pulse2SI1_Gain + UnTrans_basicLQR_B.psi;

  /* Saturate: '<S12>/Saturation' incorporates:
   *  Constant: '<S12>/Constant'
   */
  if (UnTrans_basicLQR_P.Ts > UnTrans_basicLQR_P.Saturation_UpperSat) {
    rtb_Sumi1 = UnTrans_basicLQR_P.Saturation_UpperSat;
  } else if (UnTrans_basicLQR_P.Ts < UnTrans_basicLQR_P.Saturation_LowerSat) {
    rtb_Sumi1 = UnTrans_basicLQR_P.Saturation_LowerSat;
  } else {
    rtb_Sumi1 = UnTrans_basicLQR_P.Ts;
  }

  /* DiscreteFilter: '<S6>/Discrete Filter' incorporates:
   *  Memory: '<S12>/Memory2'
   *  Product: '<S12>/Divide'
   *  Saturate: '<S12>/Saturation'
   *  Sum: '<S12>/Derivative'
   */
  rtb_Fcn1 = ((UnTrans_basicLQR_B.theta -
               UnTrans_basicLQR_DW.Memory2_PreviousInput) / rtb_Sumi1 -
              UnTrans_basicLQR_P.DiscreteFilter_DenCoef[1] *
              UnTrans_basicLQR_DW.DiscreteFilter_states) /
    UnTrans_basicLQR_P.DiscreteFilter_DenCoef[0];

  /* DiscreteFilter: '<S6>/Discrete Filter' */
  UnTrans_basicLQR_B.DiscreteFilter = UnTrans_basicLQR_P.DiscreteFilter_NumCoef
    [0] * rtb_Fcn1 + UnTrans_basicLQR_P.DiscreteFilter_NumCoef[1] *
    UnTrans_basicLQR_DW.DiscreteFilter_states;

  /* Gain: '<S6>/Pulse2SI2' incorporates:
   *  Sum: '<S6>/Add1'
   */
  UnTrans_basicLQR_B.fi = (rtb_UnitDelay1 - dt) *
    UnTrans_basicLQR_P.Pulse2SI2_Gain;

  /* Saturate: '<S13>/Saturation' incorporates:
   *  Constant: '<S13>/Constant'
   */
  if (UnTrans_basicLQR_P.Ts > UnTrans_basicLQR_P.Saturation_UpperSat_i) {
    rtb_Sumi1 = UnTrans_basicLQR_P.Saturation_UpperSat_i;
  } else if (UnTrans_basicLQR_P.Ts < UnTrans_basicLQR_P.Saturation_LowerSat_f) {
    rtb_Sumi1 = UnTrans_basicLQR_P.Saturation_LowerSat_f;
  } else {
    rtb_Sumi1 = UnTrans_basicLQR_P.Ts;
  }

  /* DiscreteFilter: '<S6>/Discrete Filter1' incorporates:
   *  Memory: '<S13>/Memory2'
   *  Product: '<S13>/Divide'
   *  Saturate: '<S13>/Saturation'
   *  Sum: '<S13>/Derivative'
   */
  dt = ((UnTrans_basicLQR_B.fi - UnTrans_basicLQR_DW.Memory2_PreviousInput_j) /
        rtb_Sumi1 - UnTrans_basicLQR_P.DiscreteFilter1_DenCoef[1] *
        UnTrans_basicLQR_DW.DiscreteFilter1_states) /
    UnTrans_basicLQR_P.DiscreteFilter1_DenCoef[0];

  /* DiscreteFilter: '<S6>/Discrete Filter1' */
  UnTrans_basicLQR_B.DiscreteFilter1 =
    UnTrans_basicLQR_P.DiscreteFilter1_NumCoef[0] * dt +
    UnTrans_basicLQR_P.DiscreteFilter1_NumCoef[1] *
    UnTrans_basicLQR_DW.DiscreteFilter1_states;

  /* Relay: '<S10>/Relay' */
  if (UnTrans_basicLQR_DW.Relay_Mode) {
    rtb_Sumi1 = UnTrans_basicLQR_P.Relay_YOn;
  } else {
    rtb_Sumi1 = UnTrans_basicLQR_P.Relay_YOff;
  }

  /* Switch: '<S1>/Switch' incorporates:
   *  Constant: '<S1>/NoCtrl'
   *  Gain: '<S11>/LQR'
   *  Relay: '<S10>/Relay'
   */
  if (rtb_Sumi1 >= UnTrans_basicLQR_P.Switch_Threshold_h) {
    rtb_Gain_a[0] = UnTrans_basicLQR_P.NoCtrl_Value[0];
    rtb_Gain_a[1] = UnTrans_basicLQR_P.NoCtrl_Value[1];
  } else {
    /* Sum: '<S1>/Sum1' incorporates:
     *  Constant: '<Root>/Step [deg]'
     *  Constant: '<Root>/Step [m]'
     *  Constant: '<Root>/Upright angle correction'
     *  Constant: '<S1>/0'
     *  Gain: '<S1>/deg2rad'
     *  Gain: '<S1>/m2rad'
     *  Gain: '<S5>/Slider Gain'
     *  Gain: '<S7>/Slider Gain'
     */
    tmp[0] = UnTrans_basicLQR_P.Translation_gain *
      UnTrans_basicLQR_P.Stepm_Value * UnTrans_basicLQR_P.m2rad_Gain;
    tmp[1] = UnTrans_basicLQR_P.u_Value;
    tmp[2] = UnTrans_basicLQR_P.Uprightanglecorrection_Value;
    tmp[3] = UnTrans_basicLQR_P.u_Value;
    tmp[4] = UnTrans_basicLQR_P.Rotation_gain * UnTrans_basicLQR_P.Stepdeg_Value
      * UnTrans_basicLQR_P.deg2rad_Gain;
    tmp[5] = UnTrans_basicLQR_P.u_Value;
    tmp_0[0] = UnTrans_basicLQR_B.theta;
    tmp_0[1] = UnTrans_basicLQR_B.DiscreteFilter;
    tmp_0[2] = UnTrans_basicLQR_B.psi;
    tmp_0[3] = UnTrans_basicLQR_B.gyroOut;
    tmp_0[4] = UnTrans_basicLQR_B.fi;
    tmp_0[5] = UnTrans_basicLQR_B.DiscreteFilter1;
    for (i = 0; i < 6; i++) {
      tmp_1[i] = tmp[i] - tmp_0[i];
    }

    /* End of Sum: '<S1>/Sum1' */
    for (i = 0; i < 2; i++) {
      rtb_Gain_a[i] = 0.0;
      for (i_0 = 0; i_0 < 6; i_0++) {
        rtb_Gain_a[i] += UnTrans_basicLQR_P.LQR_Gain[(i_0 << 1) + i] * tmp_1[i_0];
      }
    }
  }

  /* End of Switch: '<S1>/Switch' */

  /* Gain: '<Root>/Gain' */
  rtb_Gain_a[0] *= UnTrans_basicLQR_P.Gain_Gain;
  rtb_Gain_a[1] *= UnTrans_basicLQR_P.Gain_Gain;

  /* Saturate: '<Root>/Saturation1' */
  if (rtb_Gain_a[0] > UnTrans_basicLQR_P.Saturation1_UpperSat) {
    rtb_UnitConversion = UnTrans_basicLQR_P.Saturation1_UpperSat;
  } else if (rtb_Gain_a[0] < UnTrans_basicLQR_P.Saturation1_LowerSat) {
    rtb_UnitConversion = UnTrans_basicLQR_P.Saturation1_LowerSat;
  } else {
    rtb_UnitConversion = rtb_Gain_a[0];
  }

  /* End of Saturate: '<Root>/Saturation1' */

  /* Signum: '<S3>/Sign' */
  if (rtIsNaN(rtb_UnitConversion)) {
    rtb_Sumi1 = (rtNaN);
  } else if (rtb_UnitConversion < 0.0) {
    rtb_Sumi1 = -1.0;
  } else {
    rtb_Sumi1 = (rtb_UnitConversion > 0.0);
  }

  /* Sum: '<S3>/Sum' incorporates:
   *  Gain: '<S3>/Gain'
   *  Gain: '<S3>/Gain1'
   *  Signum: '<S3>/Sign'
   */
  UnTrans_basicLQR_B.Sum = UnTrans_basicLQR_P.CoulombViscousFriction2_offset *
    rtb_Sumi1 + UnTrans_basicLQR_P.CoulombViscousFriction2_gain *
    rtb_UnitConversion;

  /* Saturate: '<Root>/Saturation3' */
  if (rtb_Gain_a[1] > UnTrans_basicLQR_P.Saturation3_UpperSat) {
    rtb_UnitConversion = UnTrans_basicLQR_P.Saturation3_UpperSat;
  } else if (rtb_Gain_a[1] < UnTrans_basicLQR_P.Saturation3_LowerSat) {
    rtb_UnitConversion = UnTrans_basicLQR_P.Saturation3_LowerSat;
  } else {
    rtb_UnitConversion = rtb_Gain_a[1];
  }

  /* End of Saturate: '<Root>/Saturation3' */

  /* Signum: '<S2>/Sign' */
  if (rtIsNaN(rtb_UnitConversion)) {
    rtb_Sumi1 = (rtNaN);
  } else if (rtb_UnitConversion < 0.0) {
    rtb_Sumi1 = -1.0;
  } else {
    rtb_Sumi1 = (rtb_UnitConversion > 0.0);
  }

  /* Sum: '<S2>/Sum' incorporates:
   *  Gain: '<S2>/Gain'
   *  Gain: '<S2>/Gain1'
   *  Signum: '<S2>/Sign'
   */
  UnTrans_basicLQR_B.Sum_i = UnTrans_basicLQR_P.CoulombViscousFriction1_offset *
    rtb_Sumi1 + UnTrans_basicLQR_P.CoulombViscousFriction1_gain *
    rtb_UnitConversion;

  /* Scope: '<Root>/Ctrl' */
  {
    StructLogVar *svar = (StructLogVar *)
      UnTrans_basicLQR_DW.Ctrl_PWORK.LoggedData[0];
    LogVar *var = svar->signals.values;

    /* time */
    {
      double locTime = UnTrans_basicLQR_M->Timing.t[0]
        ;
      rt_UpdateLogVar((LogVar *)svar->time, &locTime, 0);
    }

    /* signals */
    {
      real_T up0[1];
      up0[0] = UnTrans_basicLQR_B.Sum;
      rt_UpdateLogVar((LogVar *)var, up0, 0);
      var = var->next;
    }

    {
      real_T up1[1];
      up1[0] = UnTrans_basicLQR_B.Sum_i;
      rt_UpdateLogVar((LogVar *)var, up1, 0);
    }
  }

  /* Scope: '<Root>/State' */
  {
    StructLogVar *svar = (StructLogVar *)
      UnTrans_basicLQR_DW.State_PWORK.LoggedData[0];
    LogVar *var = svar->signals.values;

    /* time */
    {
      double locTime = UnTrans_basicLQR_M->Timing.t[0]
        ;
      rt_UpdateLogVar((LogVar *)svar->time, &locTime, 0);
    }

    /* signals */
    {
      real_T up0[1];
      up0[0] = UnTrans_basicLQR_B.theta;
      rt_UpdateLogVar((LogVar *)var, up0, 0);
      var = var->next;
    }

    {
      real_T up1[1];
      up1[0] = UnTrans_basicLQR_B.DiscreteFilter;
      rt_UpdateLogVar((LogVar *)var, up1, 0);
      var = var->next;
    }

    {
      real_T up2[1];
      up2[0] = UnTrans_basicLQR_B.psi;
      rt_UpdateLogVar((LogVar *)var, up2, 0);
      var = var->next;
    }

    {
      real_T up3[1];
      up3[0] = UnTrans_basicLQR_B.gyroOut;
      rt_UpdateLogVar((LogVar *)var, up3, 0);
      var = var->next;
    }

    {
      real_T up4[1];
      up4[0] = UnTrans_basicLQR_B.fi;
      rt_UpdateLogVar((LogVar *)var, up4, 0);
      var = var->next;
    }

    {
      real_T up5[1];
      up5[0] = UnTrans_basicLQR_B.DiscreteFilter1;
      rt_UpdateLogVar((LogVar *)var, up5, 0);
    }
  }

  /* ManualSwitch: '<S4>/Reset Encoders Watchdog and Integrals' incorporates:
   *  DigitalClock: '<S4>/Digital Clock'
   *  Switch: '<S4>/Switch'
   */
  if (UnTrans_basicLQR_P.ResetEncodersWatchdogandIntegra == 1) {
    /* ManualSwitch: '<S4>/Reset Encoders Watchdog and Integrals' incorporates:
     *  Constant: '<S4>/Reset'
     */
    UnTrans_basicLQR_B.ResetEncodersWatchdogandIntegra =
      UnTrans_basicLQR_P.Reset_Value;
  } else if (UnTrans_basicLQR_M->Timing.t[0] >
             UnTrans_basicLQR_P.Switch_Threshold) {
    /* Switch: '<S4>/Switch' incorporates:
     *  Constant: '<S4>/Normal'
     *  ManualSwitch: '<S4>/Reset Encoders Watchdog and Integrals'
     */
    UnTrans_basicLQR_B.ResetEncodersWatchdogandIntegra =
      UnTrans_basicLQR_P.Normal_Value;
  } else {
    /* ManualSwitch: '<S4>/Reset Encoders Watchdog and Integrals' incorporates:
     *  Constant: '<S4>/Reset'
     */
    UnTrans_basicLQR_B.ResetEncodersWatchdogandIntegra =
      UnTrans_basicLQR_P.Reset_Value;
  }

  /* End of ManualSwitch: '<S4>/Reset Encoders Watchdog and Integrals' */

  /* Constant: '<S8>/Constant1' */
  UnTrans_basicLQR_B.PWMperiodsec = UnTrans_basicLQR_P.Constant1_Value;

  /* Gain: '<S8>/Gain' */
  UnTrans_basicLQR_B.Gain = UnTrans_basicLQR_P.Gain_Gain_n *
    UnTrans_basicLQR_B.Sum_i;

  /* Gain: '<S8>/Gain2' */
  UnTrans_basicLQR_B.Gain2 = UnTrans_basicLQR_P.Gain2_Gain *
    UnTrans_basicLQR_B.Sum;

  /* DiscreteIntegrator: '<S9>/Sumi' */
  UnTrans_basicLQR_B.Sumi = UnTrans_basicLQR_DW.Sumi_DSTATE;

  /* S-Function (microZedSwitch): '<S15>/S-Function' */

  /* Level2 S-Function Block: '<S15>/S-Function' (microZedSwitch) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[1];
    sfcnOutputs(rts,0);
  }

  /* S-Function (microZedLED): '<S15>/S-Function1' */

  /* Level2 S-Function Block: '<S15>/S-Function1' (microZedLED) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[2];
    sfcnOutputs(rts,0);
  }

  /* DiscreteIntegrator: '<S15>/Sumi1' */
  rtb_Sumi1 = UnTrans_basicLQR_DW.Sumi1_DSTATE;

  /* Lookup_n-D: '<S15>/1-D Lookup Table' incorporates:
   *  DiscreteIntegrator: '<S15>/Sumi1'
   *  Fcn: '<S15>/Fcn1'
   */
  UnTrans_basicLQR_B.uDLookupTable = look1_binlxpw(rt_remd_snf
    (UnTrans_basicLQR_DW.Sumi1_DSTATE, 1.0),
    UnTrans_basicLQR_P.uDLookupTable_bp01Data,
    UnTrans_basicLQR_P.uDLookupTable_tableData, 7U);

  /* Sum: '<S15>/Add2' incorporates:
   *  DiscreteIntegrator: '<S15>/Sumi1'
   *  Memory: '<S15>/Memory1'
   */
  rtb_Add2 = UnTrans_basicLQR_DW.Sumi1_DSTATE -
    UnTrans_basicLQR_DW.Memory1_PreviousInput;

  /* S-Function (microZedTimer): '<S15>/S-Function2' */

  /* Level2 S-Function Block: '<S15>/S-Function2' (microZedTimer) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[3];
    sfcnOutputs(rts,0);
  }

  /* Gain: '<S15>/Mu1' incorporates:
   *  Memory: '<S15>/Memory'
   *  Sum: '<S15>/Add1'
   */
  rtb_UnitDelay1 = (UnTrans_basicLQR_B.SFunction2 -
                    UnTrans_basicLQR_DW.Memory_PreviousInput) *
    UnTrans_basicLQR_P.Mu1_Gain;

  /* Outputs for Enabled SubSystem: '<S15>/Min//Max Period' incorporates:
   *  EnablePort: '<S16>/Enable'
   */
  /* Fcn: '<S15>/Fcn' incorporates:
   *  DigitalClock: '<S15>/Digital Clock'
   */
  if (UnTrans_basicLQR_M->Timing.t[0] > 0.1) {
    /* Gain: '<S16>/Mu2' */
    rtb_UnitConversion = UnTrans_basicLQR_P.Mu2_Gain * rtb_UnitDelay1;

    /* MinMax: '<S16>/MinMax' incorporates:
     *  UnitDelay: '<S16>/Unit Delay'
     */
    if ((rtb_UnitConversion <= UnTrans_basicLQR_DW.UnitDelay_DSTATE) || rtIsNaN
        (UnTrans_basicLQR_DW.UnitDelay_DSTATE)) {
      /* MinMax: '<S16>/MinMax' */
      UnTrans_basicLQR_B.MinMax = rtb_UnitConversion;
    } else {
      /* MinMax: '<S16>/MinMax' */
      UnTrans_basicLQR_B.MinMax = UnTrans_basicLQR_DW.UnitDelay_DSTATE;
    }

    /* End of MinMax: '<S16>/MinMax' */

    /* MinMax: '<S16>/MinMax1' incorporates:
     *  UnitDelay: '<S16>/Unit Delay1'
     */
    if ((rtb_UnitConversion >= UnTrans_basicLQR_DW.UnitDelay1_DSTATE) || rtIsNaN
        (UnTrans_basicLQR_DW.UnitDelay1_DSTATE)) {
      /* MinMax: '<S16>/MinMax1' */
      UnTrans_basicLQR_B.MinMax1 = rtb_UnitConversion;
    } else {
      /* MinMax: '<S16>/MinMax1' */
      UnTrans_basicLQR_B.MinMax1 = UnTrans_basicLQR_DW.UnitDelay1_DSTATE;
    }

    /* End of MinMax: '<S16>/MinMax1' */

    /* Update for UnitDelay: '<S16>/Unit Delay' */
    UnTrans_basicLQR_DW.UnitDelay_DSTATE = UnTrans_basicLQR_B.MinMax;

    /* Update for UnitDelay: '<S16>/Unit Delay1' */
    UnTrans_basicLQR_DW.UnitDelay1_DSTATE = UnTrans_basicLQR_B.MinMax1;
    srUpdateBC(UnTrans_basicLQR_DW.MinMaxPeriod_SubsysRanBC);
  }

  /* End of Fcn: '<S15>/Fcn' */
  /* End of Outputs for SubSystem: '<S15>/Min//Max Period' */

  /* Gain: '<S15>/Mu2' */
  rtb_Mu2 = UnTrans_basicLQR_P.Mu2_Gain_p * rtb_UnitDelay1;

  /* Constant: '<S15>/Constant2' */
  UnTrans_basicLQR_B.Constant2 = UnTrans_basicLQR_P.Constant2_Value;

  /* S-Function (microZedStats): '<S15>/S-Function3' */

  /* Level2 S-Function Block: '<S15>/S-Function3' (microZedStats) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[4];
    sfcnOutputs(rts,0);
  }

  /* Update for Memory: '<S12>/Memory2' */
  UnTrans_basicLQR_DW.Memory2_PreviousInput = UnTrans_basicLQR_B.theta;

  /* Update for DiscreteFilter: '<S6>/Discrete Filter' */
  UnTrans_basicLQR_DW.DiscreteFilter_states = rtb_Fcn1;

  /* Update for Memory: '<S13>/Memory2' */
  UnTrans_basicLQR_DW.Memory2_PreviousInput_j = UnTrans_basicLQR_B.fi;

  /* Update for DiscreteFilter: '<S6>/Discrete Filter1' */
  UnTrans_basicLQR_DW.DiscreteFilter1_states = dt;

  /* Update for DiscreteIntegrator: '<S9>/Sumi' incorporates:
   *  Constant: '<S9>/Constant4'
   */
  UnTrans_basicLQR_DW.Sumi_DSTATE += UnTrans_basicLQR_P.Sumi_gainval *
    UnTrans_basicLQR_P.Constant4_Value;

  /* Update for DiscreteIntegrator: '<S15>/Sumi1' incorporates:
   *  Constant: '<S15>/Constant1'
   */
  UnTrans_basicLQR_DW.Sumi1_DSTATE += UnTrans_basicLQR_P.Sumi1_gainval *
    UnTrans_basicLQR_P.Constant1_Value_f;

  /* Update for Memory: '<S15>/Memory1' */
  UnTrans_basicLQR_DW.Memory1_PreviousInput = rtb_Sumi1;

  /* Update for Memory: '<S15>/Memory' */
  UnTrans_basicLQR_DW.Memory_PreviousInput = UnTrans_basicLQR_B.SFunction2;

  /* Matfile logging */
  rt_UpdateTXYLogVars(UnTrans_basicLQR_M->rtwLogInfo,
                      (UnTrans_basicLQR_M->Timing.t));

  /* External mode */
  rtExtModeUploadCheckTrigger(1);

  {                                    /* Sample time: [0.01s, 0.0s] */
    rtExtModeUpload(0, (real_T)UnTrans_basicLQR_M->Timing.t[0]);
  }

  /* signal main to stop simulation */
  {                                    /* Sample time: [0.01s, 0.0s] */
    if ((rtmGetTFinal(UnTrans_basicLQR_M)!=-1) &&
        !((rtmGetTFinal(UnTrans_basicLQR_M)-UnTrans_basicLQR_M->Timing.t[0]) >
          UnTrans_basicLQR_M->Timing.t[0] * (DBL_EPSILON))) {
      rtmSetErrorStatus(UnTrans_basicLQR_M, "Simulation finished");
    }

    if (rtmGetStopRequested(UnTrans_basicLQR_M)) {
      rtmSetErrorStatus(UnTrans_basicLQR_M, "Simulation finished");
    }
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
  if (!(++UnTrans_basicLQR_M->Timing.clockTick0)) {
    ++UnTrans_basicLQR_M->Timing.clockTickH0;
  }

  UnTrans_basicLQR_M->Timing.t[0] = UnTrans_basicLQR_M->Timing.clockTick0 *
    UnTrans_basicLQR_M->Timing.stepSize0 +
    UnTrans_basicLQR_M->Timing.clockTickH0 *
    UnTrans_basicLQR_M->Timing.stepSize0 * 4294967296.0;
}

/* Model initialize function */
void UnTrans_basicLQR_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)UnTrans_basicLQR_M, 0,
                sizeof(RT_MODEL_UnTrans_basicLQR_T));
  rtsiSetSolverName(&UnTrans_basicLQR_M->solverInfo,"FixedStepDiscrete");
  UnTrans_basicLQR_M->solverInfoPtr = (&UnTrans_basicLQR_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap = UnTrans_basicLQR_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;

    /* polyspace +2 MISRA2012:D4.1 [Justified:Low] "UnTrans_basicLQR_M points to
       static memory which is guaranteed to be non-NULL" */
    UnTrans_basicLQR_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    UnTrans_basicLQR_M->Timing.sampleTimes =
      (&UnTrans_basicLQR_M->Timing.sampleTimesArray[0]);
    UnTrans_basicLQR_M->Timing.offsetTimes =
      (&UnTrans_basicLQR_M->Timing.offsetTimesArray[0]);

    /* task periods */
    UnTrans_basicLQR_M->Timing.sampleTimes[0] = (0.01);

    /* task offsets */
    UnTrans_basicLQR_M->Timing.offsetTimes[0] = (0.0);
  }

  rtmSetTPtr(UnTrans_basicLQR_M, &UnTrans_basicLQR_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = UnTrans_basicLQR_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    UnTrans_basicLQR_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(UnTrans_basicLQR_M, 60.0);
  UnTrans_basicLQR_M->Timing.stepSize0 = 0.01;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = (NULL);
    UnTrans_basicLQR_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(UnTrans_basicLQR_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(UnTrans_basicLQR_M->rtwLogInfo, (NULL));
    rtliSetLogT(UnTrans_basicLQR_M->rtwLogInfo, "");
    rtliSetLogX(UnTrans_basicLQR_M->rtwLogInfo, "");
    rtliSetLogXFinal(UnTrans_basicLQR_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(UnTrans_basicLQR_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(UnTrans_basicLQR_M->rtwLogInfo, 0);
    rtliSetLogMaxRows(UnTrans_basicLQR_M->rtwLogInfo, 0);
    rtliSetLogDecimation(UnTrans_basicLQR_M->rtwLogInfo, 1);
    rtliSetLogY(UnTrans_basicLQR_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(UnTrans_basicLQR_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(UnTrans_basicLQR_M->rtwLogInfo, (NULL));
  }

  /* External mode info */
  UnTrans_basicLQR_M->Sizes.checksums[0] = (3106538297U);
  UnTrans_basicLQR_M->Sizes.checksums[1] = (1166197426U);
  UnTrans_basicLQR_M->Sizes.checksums[2] = (1437391692U);
  UnTrans_basicLQR_M->Sizes.checksums[3] = (2965872559U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[7];
    UnTrans_basicLQR_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = &rtAlwaysEnabled;
    systemRan[2] = &rtAlwaysEnabled;
    systemRan[3] = &rtAlwaysEnabled;
    systemRan[4] = &rtAlwaysEnabled;
    systemRan[5] = &rtAlwaysEnabled;
    systemRan[6] = (sysRanDType *)&UnTrans_basicLQR_DW.MinMaxPeriod_SubsysRanBC;
    rteiSetModelMappingInfoPtr(UnTrans_basicLQR_M->extModeInfo,
      &UnTrans_basicLQR_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(UnTrans_basicLQR_M->extModeInfo,
                        UnTrans_basicLQR_M->Sizes.checksums);
    rteiSetTPtr(UnTrans_basicLQR_M->extModeInfo, rtmGetTPtr(UnTrans_basicLQR_M));
  }

  UnTrans_basicLQR_M->solverInfoPtr = (&UnTrans_basicLQR_M->solverInfo);
  UnTrans_basicLQR_M->Timing.stepSize = (0.01);
  rtsiSetFixedStepSize(&UnTrans_basicLQR_M->solverInfo, 0.01);
  rtsiSetSolverMode(&UnTrans_basicLQR_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  {
    int32_T i;
    for (i = 0; i < 6; i++) {
      UnTrans_basicLQR_B.IMU[i] = 0.0;
    }

    for (i = 0; i < 9; i++) {
      UnTrans_basicLQR_B.SFunction3[i] = 0.0;
    }

    UnTrans_basicLQR_B.Encoders[0] = 0.0;
    UnTrans_basicLQR_B.Encoders[1] = 0.0;
    UnTrans_basicLQR_B.PWMstatus = 0.0;
    UnTrans_basicLQR_B.IMUstatus = 0.0;
    UnTrans_basicLQR_B.theta = 0.0;
    UnTrans_basicLQR_B.DiscreteFilter = 0.0;
    UnTrans_basicLQR_B.fi = 0.0;
    UnTrans_basicLQR_B.DiscreteFilter1 = 0.0;
    UnTrans_basicLQR_B.Sum = 0.0;
    UnTrans_basicLQR_B.Sum_i = 0.0;
    UnTrans_basicLQR_B.ResetEncodersWatchdogandIntegra = 0.0;
    UnTrans_basicLQR_B.PWMperiodsec = 0.0;
    UnTrans_basicLQR_B.Gain = 0.0;
    UnTrans_basicLQR_B.Gain2 = 0.0;
    UnTrans_basicLQR_B.Sumi = 0.0;
    UnTrans_basicLQR_B.SFunction = 0.0;
    UnTrans_basicLQR_B.uDLookupTable = 0.0;
    UnTrans_basicLQR_B.SFunction2 = 0.0;
    UnTrans_basicLQR_B.Constant2 = 0.0;
    UnTrans_basicLQR_B.MinMax = 0.0;
    UnTrans_basicLQR_B.MinMax1 = 0.0;
    UnTrans_basicLQR_B.psi = 0.0;
    UnTrans_basicLQR_B.gyroOut = 0.0;
  }

  /* states (dwork) */
  (void) memset((void *)&UnTrans_basicLQR_DW, 0,
                sizeof(DW_UnTrans_basicLQR_T));
  UnTrans_basicLQR_DW.DiscreteFilter_states = 0.0;
  UnTrans_basicLQR_DW.DiscreteFilter1_states = 0.0;
  UnTrans_basicLQR_DW.Sumi_DSTATE = 0.0;
  UnTrans_basicLQR_DW.Sumi1_DSTATE = 0.0;
  UnTrans_basicLQR_DW.UnitDelay_DSTATE = 0.0;
  UnTrans_basicLQR_DW.UnitDelay1_DSTATE = 0.0;
  UnTrans_basicLQR_DW.Memory2_PreviousInput = 0.0;
  UnTrans_basicLQR_DW.Memory2_PreviousInput_j = 0.0;
  UnTrans_basicLQR_DW.Memory1_PreviousInput = 0.0;
  UnTrans_basicLQR_DW.Memory_PreviousInput = 0.0;
  UnTrans_basicLQR_DW.H[0] = 0.0;
  UnTrans_basicLQR_DW.H[1] = 0.0;
  UnTrans_basicLQR_DW.H[2] = 0.0;

  {
    int32_T i;
    for (i = 0; i < 9; i++) {
      UnTrans_basicLQR_DW.Q[i] = 0.0;
    }
  }

  UnTrans_basicLQR_DW.R = 0.0;

  {
    int32_T i;
    for (i = 0; i < 9; i++) {
      UnTrans_basicLQR_DW.P[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 9; i++) {
      UnTrans_basicLQR_DW.I[i] = 0.0;
    }
  }

  UnTrans_basicLQR_DW.x_est[0] = 0.0;
  UnTrans_basicLQR_DW.x_est[1] = 0.0;
  UnTrans_basicLQR_DW.x_est[2] = 0.0;
  UnTrans_basicLQR_DW.q = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    UnTrans_basicLQR_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 23;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo = &UnTrans_basicLQR_M->NonInlinedSFcns.sfcnInfo;
    UnTrans_basicLQR_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, (&rtmGetErrorStatus(UnTrans_basicLQR_M)));
    UnTrans_basicLQR_M->Sizes.numSampTimes = (1);
    rtssSetNumRootSampTimesPtr(sfcnInfo, &UnTrans_basicLQR_M->Sizes.numSampTimes);
    UnTrans_basicLQR_M->NonInlinedSFcns.taskTimePtrs[0] = &(rtmGetTPtr
      (UnTrans_basicLQR_M)[0]);
    rtssSetTPtrPtr(sfcnInfo,UnTrans_basicLQR_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(UnTrans_basicLQR_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(UnTrans_basicLQR_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput
      (UnTrans_basicLQR_M));
    rtssSetStepSizePtr(sfcnInfo, &UnTrans_basicLQR_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested(UnTrans_basicLQR_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &UnTrans_basicLQR_M->derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo, &UnTrans_basicLQR_M->zCCacheNeedsReset);
    rtssSetContTimeOutputInconsistentWithStateAtMajorStepPtr(sfcnInfo,
      &UnTrans_basicLQR_M->CTOutputIncnstWithState);
    rtssSetSampleHitsPtr(sfcnInfo, &UnTrans_basicLQR_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &UnTrans_basicLQR_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &UnTrans_basicLQR_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo, &UnTrans_basicLQR_M->solverInfoPtr);
  }

  UnTrans_basicLQR_M->Sizes.numSFcns = (5);

  /* register each child */
  {
    (void) memset((void *)&UnTrans_basicLQR_M->NonInlinedSFcns.childSFunctions[0],
                  0,
                  5*sizeof(SimStruct));
    UnTrans_basicLQR_M->childSfunctions =
      (&UnTrans_basicLQR_M->NonInlinedSFcns.childSFunctionPtrs[0]);

    {
      int_T i;
      for (i = 0; i < 5; i++) {
        UnTrans_basicLQR_M->childSfunctions[i] =
          (&UnTrans_basicLQR_M->NonInlinedSFcns.childSFunctions[i]);
      }
    }

    /* Level2 S-Function Block: UnTrans_basicLQR/<S8>/S-Function3 (microZedUnTrans) */
    {
      SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      {
        ssSetBlkInfo2Ptr(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.blkInfo2[0]);
      }

      _ssSetBlkInfo2PortInfo2Ptr(rts,
        &UnTrans_basicLQR_M->NonInlinedSFcns.inputOutputPortInfo2[0]);

      /* Set up the mdlInfo pointer */
      ssSetRTWSfcnInfo(rts, UnTrans_basicLQR_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods2[0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods3[0]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods4[0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.statesInfo2[0]);
        ssSetPeriodicStatesInfo(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.periodicStatesInfo[0]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 3);
        ssSetPortInfoForInputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);
        ssSetPortInfoForInputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);
        _ssSetPortInfo2ForInputUnits(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.inputPortUnits[0]);
        ssSetInputPortUnit(rts, 0, 0);
        ssSetInputPortUnit(rts, 1, 0);
        ssSetInputPortUnit(rts, 2, 0);
        _ssSetPortInfo2ForInputCoSimAttribute(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.inputPortCoSimAttribute[0]);
        ssSetInputPortIsContinuousQuantity(rts, 0, 0);
        ssSetInputPortIsContinuousQuantity(rts, 1, 0);
        ssSetInputPortIsContinuousQuantity(rts, 2, 0);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.UPtrs0;
          sfcnUPtrs[0] = &UnTrans_basicLQR_B.ResetEncodersWatchdogandIntegra;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidthAsInt(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.UPtrs1;
          sfcnUPtrs[0] = &UnTrans_basicLQR_B.PWMperiodsec;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidthAsInt(rts, 1, 1);
        }

        /* port 2 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.UPtrs2;
          sfcnUPtrs[0] = &UnTrans_basicLQR_B.Gain2;
          sfcnUPtrs[1] = &UnTrans_basicLQR_B.Gain;
          ssSetInputPortSignalPtrs(rts, 2, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 2, 1);
          ssSetInputPortWidthAsInt(rts, 2, 2);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.outputPortInfo[0]);
        ssSetPortInfoForOutputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 4);
        _ssSetPortInfo2ForOutputUnits(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.outputPortUnits[0]);
        ssSetOutputPortUnit(rts, 0, 0);
        ssSetOutputPortUnit(rts, 1, 0);
        ssSetOutputPortUnit(rts, 2, 0);
        ssSetOutputPortUnit(rts, 3, 0);
        _ssSetPortInfo2ForOutputCoSimAttribute(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.outputPortCoSimAttribute[0]);
        ssSetOutputPortIsContinuousQuantity(rts, 0, 0);
        ssSetOutputPortIsContinuousQuantity(rts, 1, 0);
        ssSetOutputPortIsContinuousQuantity(rts, 2, 0);
        ssSetOutputPortIsContinuousQuantity(rts, 3, 0);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidthAsInt(rts, 0, 2);
          ssSetOutputPortSignal(rts, 0, ((real_T *) UnTrans_basicLQR_B.Encoders));
        }

        /* port 1 */
        {
          _ssSetOutputPortNumDimensions(rts, 1, 1);
          ssSetOutputPortWidthAsInt(rts, 1, 6);
          ssSetOutputPortSignal(rts, 1, ((real_T *) UnTrans_basicLQR_B.IMU));
        }

        /* port 2 */
        {
          _ssSetOutputPortNumDimensions(rts, 2, 1);
          ssSetOutputPortWidthAsInt(rts, 2, 1);
          ssSetOutputPortSignal(rts, 2, ((real_T *)
            &UnTrans_basicLQR_B.PWMstatus));
        }

        /* port 3 */
        {
          _ssSetOutputPortNumDimensions(rts, 3, 1);
          ssSetOutputPortWidthAsInt(rts, 3, 1);
          ssSetOutputPortSignal(rts, 3, ((real_T *)
            &UnTrans_basicLQR_B.IMUstatus));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function3");
      ssSetPath(rts, "UnTrans_basicLQR/UnTrans/S-Function3");
      ssSetRTModel(rts,UnTrans_basicLQR_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn0.params;
        ssSetSFcnParamsCount(rts, 1);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)UnTrans_basicLQR_P.SFunction3_P1_Size);
      }

      /* registration */
      microZedUnTrans(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.01);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCsAsInt(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetInputPortConnected(rts, 2, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 2, 0);
      _ssSetOutputPortConnected(rts, 3, 0);
      _ssSetOutputPortBeingMerged(rts, 0, 0);
      _ssSetOutputPortBeingMerged(rts, 1, 0);
      _ssSetOutputPortBeingMerged(rts, 2, 0);
      _ssSetOutputPortBeingMerged(rts, 3, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
      ssSetInputPortBufferDstPort(rts, 2, -1);
    }

    /* Level2 S-Function Block: UnTrans_basicLQR/<S15>/S-Function (microZedSwitch) */
    {
      SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[1];

      /* timing info */
      time_T *sfcnPeriod = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn1.sfcnPeriod;
      time_T *sfcnOffset = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn1.sfcnOffset;
      int_T *sfcnTsMap = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn1.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      {
        ssSetBlkInfo2Ptr(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.blkInfo2[1]);
      }

      _ssSetBlkInfo2PortInfo2Ptr(rts,
        &UnTrans_basicLQR_M->NonInlinedSFcns.inputOutputPortInfo2[1]);

      /* Set up the mdlInfo pointer */
      ssSetRTWSfcnInfo(rts, UnTrans_basicLQR_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods2[1]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods3[1]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods4[1]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.statesInfo2[1]);
        ssSetPeriodicStatesInfo(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.periodicStatesInfo[1]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn1.outputPortInfo[0]);
        ssSetPortInfoForOutputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn1.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);
        _ssSetPortInfo2ForOutputUnits(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn1.outputPortUnits[0]);
        ssSetOutputPortUnit(rts, 0, 0);
        _ssSetPortInfo2ForOutputCoSimAttribute(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn1.outputPortCoSimAttribute[0]);
        ssSetOutputPortIsContinuousQuantity(rts, 0, 0);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidthAsInt(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &UnTrans_basicLQR_B.SFunction));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts, "UnTrans_basicLQR/microZed /MicroZed/S-Function");
      ssSetRTModel(rts,UnTrans_basicLQR_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn1.params;
        ssSetSFcnParamsCount(rts, 1);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)UnTrans_basicLQR_P.SFunction3_P1_Size);
      }

      /* registration */
      microZedSwitch(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.01);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCsAsInt(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
    }

    /* Level2 S-Function Block: UnTrans_basicLQR/<S15>/S-Function1 (microZedLED) */
    {
      SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[2];

      /* timing info */
      time_T *sfcnPeriod = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn2.sfcnPeriod;
      time_T *sfcnOffset = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn2.sfcnOffset;
      int_T *sfcnTsMap = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn2.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      {
        ssSetBlkInfo2Ptr(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.blkInfo2[2]);
      }

      _ssSetBlkInfo2PortInfo2Ptr(rts,
        &UnTrans_basicLQR_M->NonInlinedSFcns.inputOutputPortInfo2[2]);

      /* Set up the mdlInfo pointer */
      ssSetRTWSfcnInfo(rts, UnTrans_basicLQR_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods2[2]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods3[2]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods4[2]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.statesInfo2[2]);
        ssSetPeriodicStatesInfo(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.periodicStatesInfo[2]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn2.inputPortInfo[0]);
        ssSetPortInfoForInputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn2.inputPortInfo[0]);
        _ssSetPortInfo2ForInputUnits(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn2.inputPortUnits[0]);
        ssSetInputPortUnit(rts, 0, 0);
        _ssSetPortInfo2ForInputCoSimAttribute(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn2.inputPortCoSimAttribute[0]);
        ssSetInputPortIsContinuousQuantity(rts, 0, 0);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn2.UPtrs0;
          sfcnUPtrs[0] = &UnTrans_basicLQR_B.uDLookupTable;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidthAsInt(rts, 0, 1);
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function1");
      ssSetPath(rts, "UnTrans_basicLQR/microZed /MicroZed/S-Function1");
      ssSetRTModel(rts,UnTrans_basicLQR_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn2.params;
        ssSetSFcnParamsCount(rts, 1);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)UnTrans_basicLQR_P.SFunction3_P1_Size);
      }

      /* registration */
      microZedLED(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.01);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCsAsInt(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }

    /* Level2 S-Function Block: UnTrans_basicLQR/<S15>/S-Function2 (microZedTimer) */
    {
      SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[3];

      /* timing info */
      time_T *sfcnPeriod = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn3.sfcnPeriod;
      time_T *sfcnOffset = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn3.sfcnOffset;
      int_T *sfcnTsMap = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn3.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      {
        ssSetBlkInfo2Ptr(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.blkInfo2[3]);
      }

      _ssSetBlkInfo2PortInfo2Ptr(rts,
        &UnTrans_basicLQR_M->NonInlinedSFcns.inputOutputPortInfo2[3]);

      /* Set up the mdlInfo pointer */
      ssSetRTWSfcnInfo(rts, UnTrans_basicLQR_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods2[3]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods3[3]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods4[3]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.statesInfo2[3]);
        ssSetPeriodicStatesInfo(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.periodicStatesInfo[3]);
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn3.outputPortInfo[0]);
        ssSetPortInfoForOutputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn3.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);
        _ssSetPortInfo2ForOutputUnits(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn3.outputPortUnits[0]);
        ssSetOutputPortUnit(rts, 0, 0);
        _ssSetPortInfo2ForOutputCoSimAttribute(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn3.outputPortCoSimAttribute[0]);
        ssSetOutputPortIsContinuousQuantity(rts, 0, 0);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidthAsInt(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &UnTrans_basicLQR_B.SFunction2));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function2");
      ssSetPath(rts, "UnTrans_basicLQR/microZed /MicroZed/S-Function2");
      ssSetRTModel(rts,UnTrans_basicLQR_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn3.params;
        ssSetSFcnParamsCount(rts, 1);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)UnTrans_basicLQR_P.SFunction3_P1_Size);
      }

      /* registration */
      microZedTimer(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.01);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCsAsInt(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
    }

    /* Level2 S-Function Block: UnTrans_basicLQR/<S15>/S-Function3 (microZedStats) */
    {
      SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[4];

      /* timing info */
      time_T *sfcnPeriod = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.sfcnPeriod;
      time_T *sfcnOffset = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.sfcnOffset;
      int_T *sfcnTsMap = UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      {
        ssSetBlkInfo2Ptr(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.blkInfo2[4]);
      }

      _ssSetBlkInfo2PortInfo2Ptr(rts,
        &UnTrans_basicLQR_M->NonInlinedSFcns.inputOutputPortInfo2[4]);

      /* Set up the mdlInfo pointer */
      ssSetRTWSfcnInfo(rts, UnTrans_basicLQR_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods2[4]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods3[4]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.methods4[4]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &UnTrans_basicLQR_M->NonInlinedSFcns.statesInfo2[4]);
        ssSetPeriodicStatesInfo(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.periodicStatesInfo[4]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.inputPortInfo[0]);
        ssSetPortInfoForInputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.inputPortInfo[0]);
        _ssSetPortInfo2ForInputUnits(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.inputPortUnits[0]);
        ssSetInputPortUnit(rts, 0, 0);
        _ssSetPortInfo2ForInputCoSimAttribute(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.inputPortCoSimAttribute[0]);
        ssSetInputPortIsContinuousQuantity(rts, 0, 0);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.UPtrs0;
          sfcnUPtrs[0] = &UnTrans_basicLQR_B.Constant2;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidthAsInt(rts, 0, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.outputPortInfo[0]);
        ssSetPortInfoForOutputs(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);
        _ssSetPortInfo2ForOutputUnits(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.outputPortUnits[0]);
        ssSetOutputPortUnit(rts, 0, 0);
        _ssSetPortInfo2ForOutputCoSimAttribute(rts,
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.outputPortCoSimAttribute[0]);
        ssSetOutputPortIsContinuousQuantity(rts, 0, 0);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidthAsInt(rts, 0, 9);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            UnTrans_basicLQR_B.SFunction3));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function3");
      ssSetPath(rts, "UnTrans_basicLQR/microZed /MicroZed/S-Function3");
      ssSetRTModel(rts,UnTrans_basicLQR_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &UnTrans_basicLQR_M->NonInlinedSFcns.Sfcn4.params;
        ssSetSFcnParamsCount(rts, 1);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)UnTrans_basicLQR_P.SFunction3_P1_Size);
      }

      /* registration */
      microZedStats(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.01);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCsAsInt(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }
  }

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(UnTrans_basicLQR_M->rtwLogInfo, 0.0,
    rtmGetTFinal(UnTrans_basicLQR_M), UnTrans_basicLQR_M->Timing.stepSize0,
    (&rtmGetErrorStatus(UnTrans_basicLQR_M)));

  /* SetupRuntimeResources for Scope: '<Root>/Ctrl' */
  {
    RTWLogSignalInfo rt_ScopeSignalInfo;
    static int_T rt_ScopeSignalWidths[] = { 1, 1 };

    static int_T rt_ScopeSignalNumDimensions[] = { 1, 1 };

    static int_T rt_ScopeSignalDimensions[] = { 1, 1 };

    static void *rt_ScopeCurrSigDims[] = { (NULL), (NULL) };

    static int_T rt_ScopeCurrSigDimsSize[] = { 4, 4 };

    static const char_T *rt_ScopeSignalLabels[] = { "Left Control",
      "Right Control" };

    static char_T rt_ScopeSignalTitles[] = "Left ControlRight Control";
    static int_T rt_ScopeSignalTitleLengths[] = { 12, 13 };

    static boolean_T rt_ScopeSignalIsVarDims[] = { 0, 0 };

    static int_T rt_ScopeSignalPlotStyles[] = { 1, 1 };

    BuiltInDTypeId dTypes[2] = { SS_DOUBLE, SS_DOUBLE };

    static char_T rt_ScopeBlockName[] = "UnTrans_basicLQR/Ctrl";
    static int_T rt_ScopeFrameData[] = { 0, 0 };

    static RTWPreprocessingFcnPtr rt_ScopeSignalLoggingPreprocessingFcnPtrs[] =
      {
      (NULL), (NULL)
    };

    rt_ScopeSignalInfo.numSignals = 2;
    rt_ScopeSignalInfo.numCols = rt_ScopeSignalWidths;
    rt_ScopeSignalInfo.numDims = rt_ScopeSignalNumDimensions;
    rt_ScopeSignalInfo.dims = rt_ScopeSignalDimensions;
    rt_ScopeSignalInfo.isVarDims = rt_ScopeSignalIsVarDims;
    rt_ScopeSignalInfo.currSigDims = rt_ScopeCurrSigDims;
    rt_ScopeSignalInfo.currSigDimsSize = rt_ScopeCurrSigDimsSize;
    rt_ScopeSignalInfo.dataTypes = dTypes;
    rt_ScopeSignalInfo.complexSignals = (NULL);
    rt_ScopeSignalInfo.frameData = rt_ScopeFrameData;
    rt_ScopeSignalInfo.preprocessingPtrs =
      rt_ScopeSignalLoggingPreprocessingFcnPtrs;
    rt_ScopeSignalInfo.labels.cptr = rt_ScopeSignalLabels;
    rt_ScopeSignalInfo.titles = rt_ScopeSignalTitles;
    rt_ScopeSignalInfo.titleLengths = rt_ScopeSignalTitleLengths;
    rt_ScopeSignalInfo.plotStyles = rt_ScopeSignalPlotStyles;
    rt_ScopeSignalInfo.blockNames.cptr = (NULL);
    rt_ScopeSignalInfo.stateNames.cptr = (NULL);
    rt_ScopeSignalInfo.crossMdlRef = (NULL);
    rt_ScopeSignalInfo.dataTypeConvert = (NULL);
    UnTrans_basicLQR_DW.Ctrl_PWORK.LoggedData[0] = rt_CreateStructLogVar(
      UnTrans_basicLQR_M->rtwLogInfo,
      0.0,
      rtmGetTFinal(UnTrans_basicLQR_M),
      UnTrans_basicLQR_M->Timing.stepSize0,
      (&rtmGetErrorStatus(UnTrans_basicLQR_M)),
      "CONTROL",
      1,
      3000,
      1,
      0.01,
      &rt_ScopeSignalInfo,
      rt_ScopeBlockName);
    if (UnTrans_basicLQR_DW.Ctrl_PWORK.LoggedData[0] == (NULL))
      return;
  }

  /* SetupRuntimeResources for Scope: '<Root>/State' */
  {
    RTWLogSignalInfo rt_ScopeSignalInfo;
    static int_T rt_ScopeSignalWidths[] = { 1, 1, 1, 1, 1, 1 };

    static int_T rt_ScopeSignalNumDimensions[] = { 1, 1, 1, 1, 1, 1 };

    static int_T rt_ScopeSignalDimensions[] = { 1, 1, 1, 1, 1, 1 };

    static void *rt_ScopeCurrSigDims[] = { (NULL), (NULL), (NULL), (NULL), (NULL),
      (NULL) };

    static int_T rt_ScopeCurrSigDimsSize[] = { 4, 4, 4, 4, 4, 4 };

    static const char_T *rt_ScopeSignalLabels[] = { "theta",
      "thetadot",
      "psi",
      "psidot",
      "phi",
      "phidot" };

    static char_T rt_ScopeSignalTitles[] = "thetathetadotpsipsidotphiphidot";
    static int_T rt_ScopeSignalTitleLengths[] = { 5, 8, 3, 6, 3, 6 };

    static boolean_T rt_ScopeSignalIsVarDims[] = { 0, 0, 0, 0, 0, 0 };

    static int_T rt_ScopeSignalPlotStyles[] = { 1, 1, 1, 1, 1, 1 };

    BuiltInDTypeId dTypes[6] = { SS_DOUBLE, SS_DOUBLE, SS_DOUBLE, SS_DOUBLE,
      SS_DOUBLE, SS_DOUBLE };

    static char_T rt_ScopeBlockName[] = "UnTrans_basicLQR/State";
    static int_T rt_ScopeFrameData[] = { 0, 0, 0, 0, 0, 0 };

    static RTWPreprocessingFcnPtr rt_ScopeSignalLoggingPreprocessingFcnPtrs[] =
      {
      (NULL), (NULL), (NULL), (NULL), (NULL), (NULL)
    };

    rt_ScopeSignalInfo.numSignals = 6;
    rt_ScopeSignalInfo.numCols = rt_ScopeSignalWidths;
    rt_ScopeSignalInfo.numDims = rt_ScopeSignalNumDimensions;
    rt_ScopeSignalInfo.dims = rt_ScopeSignalDimensions;
    rt_ScopeSignalInfo.isVarDims = rt_ScopeSignalIsVarDims;
    rt_ScopeSignalInfo.currSigDims = rt_ScopeCurrSigDims;
    rt_ScopeSignalInfo.currSigDimsSize = rt_ScopeCurrSigDimsSize;
    rt_ScopeSignalInfo.dataTypes = dTypes;
    rt_ScopeSignalInfo.complexSignals = (NULL);
    rt_ScopeSignalInfo.frameData = rt_ScopeFrameData;
    rt_ScopeSignalInfo.preprocessingPtrs =
      rt_ScopeSignalLoggingPreprocessingFcnPtrs;
    rt_ScopeSignalInfo.labels.cptr = rt_ScopeSignalLabels;
    rt_ScopeSignalInfo.titles = rt_ScopeSignalTitles;
    rt_ScopeSignalInfo.titleLengths = rt_ScopeSignalTitleLengths;
    rt_ScopeSignalInfo.plotStyles = rt_ScopeSignalPlotStyles;
    rt_ScopeSignalInfo.blockNames.cptr = (NULL);
    rt_ScopeSignalInfo.stateNames.cptr = (NULL);
    rt_ScopeSignalInfo.crossMdlRef = (NULL);
    rt_ScopeSignalInfo.dataTypeConvert = (NULL);
    UnTrans_basicLQR_DW.State_PWORK.LoggedData[0] = rt_CreateStructLogVar(
      UnTrans_basicLQR_M->rtwLogInfo,
      0.0,
      rtmGetTFinal(UnTrans_basicLQR_M),
      UnTrans_basicLQR_M->Timing.stepSize0,
      (&rtmGetErrorStatus(UnTrans_basicLQR_M)),
      "STATE",
      1,
      3000,
      1,
      0.01,
      &rt_ScopeSignalInfo,
      rt_ScopeBlockName);
    if (UnTrans_basicLQR_DW.State_PWORK.LoggedData[0] == (NULL))
      return;
  }

  /* Start for S-Function (microZedUnTrans): '<S8>/S-Function3' */
  /* Level2 S-Function Block: '<S8>/S-Function3' (microZedUnTrans) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[0];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for Constant: '<S8>/Constant1' */
  UnTrans_basicLQR_B.PWMperiodsec = UnTrans_basicLQR_P.Constant1_Value;

  /* Start for S-Function (microZedSwitch): '<S15>/S-Function' */
  /* Level2 S-Function Block: '<S15>/S-Function' (microZedSwitch) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[1];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for S-Function (microZedLED): '<S15>/S-Function1' */
  /* Level2 S-Function Block: '<S15>/S-Function1' (microZedLED) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[2];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for S-Function (microZedTimer): '<S15>/S-Function2' */
  /* Level2 S-Function Block: '<S15>/S-Function2' (microZedTimer) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[3];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* Start for Constant: '<S15>/Constant2' */
  UnTrans_basicLQR_B.Constant2 = UnTrans_basicLQR_P.Constant2_Value;

  /* Start for S-Function (microZedStats): '<S15>/S-Function3' */
  /* Level2 S-Function Block: '<S15>/S-Function3' (microZedStats) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[4];
    sfcnStart(rts);
    if (ssGetErrorStatus(rts) != (NULL))
      return;
  }

  /* InitializeConditions for Memory: '<S12>/Memory2' */
  UnTrans_basicLQR_DW.Memory2_PreviousInput =
    UnTrans_basicLQR_P.Memory2_InitialCondition;

  /* InitializeConditions for DiscreteFilter: '<S6>/Discrete Filter' */
  UnTrans_basicLQR_DW.DiscreteFilter_states =
    UnTrans_basicLQR_P.DiscreteFilter_InitialStates;

  /* InitializeConditions for Memory: '<S13>/Memory2' */
  UnTrans_basicLQR_DW.Memory2_PreviousInput_j =
    UnTrans_basicLQR_P.Memory2_InitialCondition_i;

  /* InitializeConditions for DiscreteFilter: '<S6>/Discrete Filter1' */
  UnTrans_basicLQR_DW.DiscreteFilter1_states =
    UnTrans_basicLQR_P.DiscreteFilter1_InitialStates;

  /* InitializeConditions for DiscreteIntegrator: '<S9>/Sumi' */
  UnTrans_basicLQR_DW.Sumi_DSTATE = UnTrans_basicLQR_P.Sumi_IC;

  /* InitializeConditions for DiscreteIntegrator: '<S15>/Sumi1' */
  UnTrans_basicLQR_DW.Sumi1_DSTATE = UnTrans_basicLQR_P.Sumi1_IC;

  /* InitializeConditions for Memory: '<S15>/Memory1' */
  UnTrans_basicLQR_DW.Memory1_PreviousInput =
    UnTrans_basicLQR_P.Memory1_InitialCondition;

  /* InitializeConditions for Memory: '<S15>/Memory' */
  UnTrans_basicLQR_DW.Memory_PreviousInput =
    UnTrans_basicLQR_P.Memory_InitialCondition;

  /* SystemInitialize for MATLAB Function: '<S6>/Kalman' */
  UnTrans_basicLQR_DW.H[0] = 1.0;
  UnTrans_basicLQR_DW.H[1] = 0.0;
  UnTrans_basicLQR_DW.H[2] = 0.0;
  UnTrans_basicLQR_DW.P_not_empty = false;
  UnTrans_basicLQR_DW.q = 4.11E-5;
  UnTrans_basicLQR_DW.R = 0.11596;

  /* SystemInitialize for Enabled SubSystem: '<S15>/Min//Max Period' */
  /* InitializeConditions for UnitDelay: '<S16>/Unit Delay' */
  UnTrans_basicLQR_DW.UnitDelay_DSTATE =
    UnTrans_basicLQR_P.UnitDelay_InitialCondition;

  /* InitializeConditions for UnitDelay: '<S16>/Unit Delay1' */
  UnTrans_basicLQR_DW.UnitDelay1_DSTATE =
    UnTrans_basicLQR_P.UnitDelay1_InitialCondition;

  /* End of SystemInitialize for SubSystem: '<S15>/Min//Max Period' */
}

/* Model terminate function */
void UnTrans_basicLQR_terminate(void)
{
  /* Terminate for S-Function (microZedUnTrans): '<S8>/S-Function3' */
  /* Level2 S-Function Block: '<S8>/S-Function3' (microZedUnTrans) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[0];
    sfcnTerminate(rts);
  }

  /* Terminate for S-Function (microZedSwitch): '<S15>/S-Function' */
  /* Level2 S-Function Block: '<S15>/S-Function' (microZedSwitch) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[1];
    sfcnTerminate(rts);
  }

  /* Terminate for S-Function (microZedLED): '<S15>/S-Function1' */
  /* Level2 S-Function Block: '<S15>/S-Function1' (microZedLED) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[2];
    sfcnTerminate(rts);
  }

  /* Terminate for S-Function (microZedTimer): '<S15>/S-Function2' */
  /* Level2 S-Function Block: '<S15>/S-Function2' (microZedTimer) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[3];
    sfcnTerminate(rts);
  }

  /* Terminate for S-Function (microZedStats): '<S15>/S-Function3' */
  /* Level2 S-Function Block: '<S15>/S-Function3' (microZedStats) */
  {
    SimStruct *rts = UnTrans_basicLQR_M->childSfunctions[4];
    sfcnTerminate(rts);
  }
}

#include <stdio.h>

/* Final time from "Simulation Parameters"   */
real_T finaltime = 60.0;

////////////////////////////////////////////////
//
//  Return compilation date and time
//
char *GetDateAndTime( void )
{
  static char AuxStr[ 128 ];
  sprintf( AuxStr, "%s %s", __DATE__, __TIME__ );
  return( AuxStr );
}

unsigned Date2Hex( void )
{
  char buff[32];
  unsigned month, day, year, d2h;
  static const char month_names[] = "JanFebMarAprMayJunJulAugSepOctNovDec";
  sscanf(__DATE__, "%s %d %d", buff, &day, &year);
  month = (strstr(month_names, buff)-month_names)/3+1;
  d2h = 0x10000000 * (year / 1000);
  year= (year % 1000);
  d2h += 0x01000000 * (year / 100);
  year= (year % 100);
  d2h += 0x00100000 * (year / 10);
  year= (year % 10);
  d2h += 0x00010000 * year;
  d2h += 0x00001000 * (month / 10);
  month= (month % 10);
  d2h += 0x00000100 * month;
  d2h += 0x00000010 * (day / 10);
  day= (day % 10);
  d2h += 0x00000001 * day;
  return( d2h );
}
