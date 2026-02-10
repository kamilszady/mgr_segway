/*
 * UnTrans_basicLQR_dt.h
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

#include "ext_types.h"

/* data type size table */
static uint_T rtDataTypeSizes[] = {
  sizeof(real_T),
  sizeof(real32_T),
  sizeof(int8_T),
  sizeof(uint8_T),
  sizeof(int16_T),
  sizeof(uint16_T),
  sizeof(int32_T),
  sizeof(uint32_T),
  sizeof(boolean_T),
  sizeof(fcn_call_T),
  sizeof(int_T),
  sizeof(pointer_T),
  sizeof(action_T),
  2*sizeof(uint32_T),
  sizeof(int32_T),
  sizeof(int64_T),
  sizeof(uint64_T),
  sizeof(uint64_T),
  sizeof(int64_T),
  sizeof(uint_T),
  sizeof(char_T),
  sizeof(uchar_T),
  sizeof(time_T)
};

/* data type name table */
static const char_T * rtDataTypeNames[] = {
  "real_T",
  "real32_T",
  "int8_T",
  "uint8_T",
  "int16_T",
  "uint16_T",
  "int32_T",
  "uint32_T",
  "boolean_T",
  "fcn_call_T",
  "int_T",
  "pointer_T",
  "action_T",
  "timer_uint32_pair_T",
  "physical_connection",
  "int64_T",
  "uint64_T",
  "uint64_T",
  "int64_T",
  "uint_T",
  "char_T",
  "uchar_T",
  "time_T"
};

/* data type transitions for block I/O structure */
static DataTypeTransition rtBTransitions[] = {
  { (char_T *)(&UnTrans_basicLQR_B.Encoders[0]), 0, 0, 38 }
  ,

  { (char_T *)(&UnTrans_basicLQR_DW.DiscreteFilter_states), 0, 0, 45 },

  { (char_T *)(&UnTrans_basicLQR_DW.Ctrl_PWORK.LoggedData[0]), 11, 0, 8 },

  { (char_T *)(&UnTrans_basicLQR_DW.MinMaxPeriod_SubsysRanBC), 2, 0, 1 },

  { (char_T *)(&UnTrans_basicLQR_DW.Relay_Mode), 8, 0, 2 }
};

/* data type transition table for block I/O structure */
static DataTypeTransitionTable rtBTransTable = {
  5U,
  rtBTransitions
};

/* data type transitions for Parameters structure */
static DataTypeTransition rtPTransitions[] = {
  { (char_T *)(&UnTrans_basicLQR_P.SFunction3_P1_Size[0]), 0, 0, 99 },

  { (char_T *)(&UnTrans_basicLQR_P.ResetEncodersWatchdogandIntegra), 3, 0, 1 }
};

/* data type transition table for Parameters structure */
static DataTypeTransitionTable rtPTransTable = {
  2U,
  rtPTransitions
};

/* [EOF] UnTrans_basicLQR_dt.h */
