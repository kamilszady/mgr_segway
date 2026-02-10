/*
 * UnTrans_basicLQR_private.h
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

#ifndef RTW_HEADER_UnTrans_basicLQR_private_h_
#define RTW_HEADER_UnTrans_basicLQR_private_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#include "UnTrans_basicLQR_types.h"
#include "UnTrans_basicLQR.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmSetTFinal
#define rtmSetTFinal(rtm, val)         ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmSetTPtr
#define rtmSetTPtr(rtm, val)           ((rtm)->Timing.t = (val))
#endif

extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern real_T rt_remd_snf(real_T u0, real_T u1);
extern real_T look1_binlxpw(real_T u0, const real_T bp0[], const real_T table[],
  uint32_T maxIndex);
extern void microZedUnTrans(SimStruct *rts);
extern void microZedSwitch(SimStruct *rts);
extern void microZedLED(SimStruct *rts);
extern void microZedTimer(SimStruct *rts);
extern void microZedStats(SimStruct *rts);

#endif                              /* RTW_HEADER_UnTrans_basicLQR_private_h_ */
