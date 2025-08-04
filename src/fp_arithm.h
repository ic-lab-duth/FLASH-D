#ifndef __FP_ARITH_H__
#define __FP_ARITH_H__

#include "ac_std_float.h"
#include "fast_float.h"

template<typename T>
inline
T m_ac_recip_pwl(const T &x) {
  T res;
  ac_math::ac_reciprocal_pwl(
    x,
    res
  );
  return T(res);
}

template<>
inline
ac::bfloat16 m_ac_recip_pwl(const ac::bfloat16 &x) {
  auto res = ac::bfloat16::zero().to_ac_std_float();
  ac_math::ac_reciprocal_pwl(
    x.to_ac_std_float(),
    res
  );
  return ac::bfloat16(res);
}

template<typename T>
inline
T m_ac_exp_pwl(const T &x) {
  T res;
  ac_math::ac_exp_pwl(
    x,
    res
  );
  return T(res);
}

template<>
inline
ac::bfloat16 m_ac_exp_pwl(const ac::bfloat16 &x) {
  auto res = ac::bfloat16::zero().to_ac_std_float();
  ac_math::ac_exp_pwl(
    x.to_ac_std_float(),
    res
  );
  return ac::bfloat16(res);
}

template<typename T>
inline
T m_mitch_ln_pwl_(const T &x) {
  static const ac_fixed<8, 1, false> 
  f = 0.693147f;
  
  static const ac_int<16, false> 
  one_bits = 0x3f80;

  ac_int<16, false> x_bits = x.data_ac_int();

  ac_fixed<16, 9, true>
  log2x = 0;
  log2x.set_slc(0, 
    ac_int<16, false>(x_bits - one_bits).template slc<16>(0)
  );

  ac_fixed<16, 9, true> lnx = f * log2x;

  return T(lnx);
}

template<typename T>
inline
T m_fma(const T A, const T B, const T C) {
  return A.template fma<AC_TRN_ZERO, false>(B, C);
}

template<typename T>
inline
T m_ln_pwl(const T &x) {
  ac_int<16,false> 
  fp_bits = x.data_ac_int();

  ac_int<8,false> exp = fp_bits.template slc<8>(7);
  ac_int<3,false> idx;
  ac_int<1,false> is_zero = !exp.or_reduce();

  switch(exp) {
    case 126: {idx = 7; break;}
    case 125: {idx = 6; break;}
    case 124: {idx = 5; break;}
    case 123: {idx = 4; break;}
    case 122: {idx = 3; break;}
    case 121: {idx = 2; break;}
    case 119: {idx = 1; break;}
    default: {idx = 0; break;}
  }

  static const T a_rom[8] = {
    T(1632.29476138034),
    T(90.0866047160712),
    T(42.3306902340338),
    T(22.5667931583203),
    T(11.0725799291238),
    T(5.53617770967981),
    T(2.78381758917953),
    T(1.37067062318228)  
  };

  static const T b_rom[8] = {
    T(-11.5292484125828),
    T(-5.50499780061298),
    T(-4.75881163683115),
    T(-4.1411898532151 ),
    T(-3.42280152639032),
    T(-2.73075124895982),
    T(-2.04266121883475),
    T(-1.33608773583612)
  };

  return m_fma(x, is_zero ? T(0.0f) : a_rom[idx], is_zero ? T(0.0f) : b_rom[idx]);
}

template<typename T>
inline
T m_sigmoid_pwl(const T &x) {
  ac_int<16,false> 
  fp_bits = x.data_ac_int();

  ac_int<1,false> sig = fp_bits[15];
  ac_int<8,false> exp = fp_bits.template slc<8>(7);
  ac_int<7,false> man = fp_bits.template slc<7>(0);

  ac_int<1, false> lt_5 = (exp == 129); 
  ac_int<1, false> lt_3 = (exp == 128) & man[6]; 
  ac_int<1, false> lt_2 = (exp == 128) & ~man[6]; 
  ac_int<1, false> lt_1 = (exp == 127);

  ac_int<3, false>
  idx = sig ? (lt_3 ? 0 : lt_2 ? 1 : lt_1 ? 2 : 0) : 
              (lt_1 ? 3 : lt_2 ? 4 : lt_3 ? 5 : lt_5 ? 6 : 7);

  static const T a_rom[8] = {
    T(0.0196664833794307 ), 
    T(0.0720769586083984 ), 
    T(0.147704383084591  ),  
    T(0.23900883040372   ),  
    T(0.14905610871022   ),  
    T(0.0662155958524261 ), 
    T(0.0356390265274004 ), 
    T(0.00226757988060892) 
  };

  static const T b_rom[8] = {
    T(0.100133753532174),
    T(0.257365179219077),
    T(0.408620028171462),
    T(0.499924475490591),
    T(0.58987719718409 ),
    T(0.755558222899678),
    T(0.847287930874755),
    T(0.980773717461921)
  };

  return m_fma(x, a_rom[idx], b_rom[idx]);
}

template<typename T, typename Tff, int N>
inline
T m_dotProd(const T A[N], const T B[N]) {
  Tff A_ff[N];
  Tff B_ff[N];

  #pragma hls_unroll yes
  BF2FF: for (int i = 0; i < N; i++) {
    A_ff[i] = from_ac2ff(A[i]);
    B_ff[i] = from_ac2ff(B[i]);
  }

  Tff dot; dot.template dotProd<N>(A_ff, B_ff);

  return from_ff2ac(dot);
}

#endif