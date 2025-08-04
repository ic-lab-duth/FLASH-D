#ifndef __DEFINES_HH__
#define __DEFINES_HH__

#include <string>
#include <ac_fixed.h>
#include <ac_std_float.h>
#include <iostream>

#include "file_io.h"
#include "logging.h"
#include "reduction.h"
#include "ac_math.h"

#define _Tc 2
#define _Bc 2
#define _d 3
#define SeqLen _Tc * _Bc
#define BF16

#ifdef FP32 
  #define fptype ac_std_float<32,8>
  #define ffptype ffloat32_t
#endif
#ifdef BF16 
  #define fptype ac::bfloat16
  #define ffptype fbfloat16_t
#endif
#ifdef FP834 
  #define fptype ac_std_float<8,4>
  #define ffptype fofp34_t
#endif
#ifdef FP852
  #define fptype ac_std_float<8,5>
  #define ffptype fofp52_t
#endif

typedef ac_int<1, false> bit_t;

struct vec_t {
  fptype data[_d];
};

struct ol2oa_ctrl {
  ac_int<1,false> in_valid;
};

#define __ONE_IEEE754(m) ac_int<1+E+M, false>(0x3f80 << (m+1-8))

#define dtype ac_fixed<_W, _I, true>

#define __AC_MIN(x) x().template set_val<AC_VAL_MIN>()
#define __AC_QNT(x) x().template set_val<AC_VAL_QUANTUM>()

#define CLOG2(x) ac::log2_ceil<x>::val
#define BITS(x)  ac::nbits<x>::val

#include "fast_float.h"

inline
fbfloat16_t
from_ac2ff(
  const ac::bfloat16 &in
){
  ac_int<16, false> in_bits = in.data_ac_int();
  return fbfloat16_t(in_bits);
}

inline
ac::bfloat16
from_ff2ac(
  const fbfloat16_t &in
){
  ac_int<16, true> in_bits = 0;

  ac::bfloat16 res;

  in_bits[15] = in.sign;
  in_bits.set_slc(7, in.exponent);
  in_bits.set_slc(0, in.mantissa);

  res.set_data(in_bits);
  return res;
}

inline
fofp34_t
from_ac2ff(
  const ac_std_float<8,4> &in
){
  ac_int<8, false> in_bits = in.data_ac_int();
  return fofp34_t(in_bits);
}

inline
ac_std_float<8,4>
from_ff2ac(
  const fofp34_t &in
){
  ac_int<8, true> in_bits = 0;

  ac_std_float<8,4> res;

  in_bits[7] = in.sign;
  in_bits.set_slc(3, in.exponent);
  in_bits.set_slc(0, in.mantissa);

  res.set_data(in_bits);
  return res;
}

inline
fofp52_t
from_ac2ff(
  const ac_std_float<8,5> &in
){
  ac_int<8, false> in_bits = in.data_ac_int();
  return fofp52_t(in_bits);
}

inline
ac_std_float<8,5>
from_ff2ac(
  const fofp52_t &in
){
  ac_int<8, true> in_bits = 0;

  ac_std_float<8,5> res;

  in_bits[7] = in.sign;
  in_bits.set_slc(2, in.exponent);
  in_bits.set_slc(0, in.mantissa);

  res.set_data(in_bits);
  return res;
}

inline
ffloat32_t
from_ac2ff(
  const ac_std_float<32,8> &in
){
  ac_int<32, false> in_bits = in.data_ac_int();
  return ffloat32_t(in_bits);
}

inline
ac_std_float<32,8>
from_ff2ac(
  const ffloat32_t &in
){
  ac_int<32, true> in_bits = 0;

  ac_std_float<32,8> res;

  in_bits[31] = in.sign;
  in_bits.set_slc(23, in.exponent);
  in_bits.set_slc(0, in.mantissa);

  res.set_data(in_bits);
  return res;
}
#endif