#ifndef __MATH_OPS_H__
#define __MATH_OPS_H__

#include "defines.h"

namespace mops{
  template< int W, int I >
  void split(ac_int<I-1, 0>     &_int,
            ac_fixed<W-I, 1, 1> &_frac,
      const ac_fixed<W, I, 1>   &in
  ){
    ac_fixed<W, I, 1> in_eff;
                      in_eff = in;

    // Compute the 2's complement. Invert and add
    // a 1 to the LSB. Because its fixed-point rep,
    // the LSB's weight is (1 >> (W-I)).
    using T = ac_fixed<W, I, 0>;
    in_eff = ac_fixed<W, I, 0> (
      ~in_eff + __AC_QNT(T)
    );

    _int.set_slc(0, in_eff.template slc<I-1>( (unsigned)(W-I) ));
    
    _frac = ac_fixed<W-I, 1, 1>(in + _int);
  }

  template< int W, int I >
  void pow2_pwl(
          ac_fixed<W-I+1, 1, 0> &a,
          ac_fixed<W-I+1, 1, 0> &b,
    const ac_fixed<W-I+1, 1, 1> &x
  ){  
    /**
     * @note If x is zero, just output 1 no need to pwl.
    */

    const ac_fixed<W-I+1, 1, 0>
    alpha_lut[8] = {
      0.663764778033316,
      0.608423508606172,
      0.557994093565377,
      0.511664863403061,
      0.469203088566355,
      0.430261502403925,
      0.394545500518956,
      0.361823424728051
    };

    const ac_fixed<W-I+1, 1, 0> 
    beta_lut[8] = {
      0.999395521280818,
      0.992477862602425,
      0.979870508842226,
      0.962497047531358,
      0.941266160113005,
      0.916927668761487,
      0.89014066734776,
      0.861508851030718
    };

    /**
      * @param index is the 3 MSBs of the fraction.
      * Here x is <1>.<W-I> --> W-I+1 bits in total
    */
    ac_int<3, false> index;
                     index.set_slc(0, x.template slc<3>(W-I+1 - 4));
                    //  index[2] = ~x[W-I];
    
    ac_int<1, false> is_zero = (x == 0);

    a = is_zero ? ac_fixed<W-I+1, 1, 0>(0.0f) : alpha_lut[index];
    b = is_zero ? ac_fixed<W-I+1, 1, 0>(1.0f) : beta_lut[index] ;
  }
  
  // I is addition integer bits + sign
  template<int W=16, int I=5>
  void
  exp_pwl(
    const ac_fixed<W, I, true> &x,
          ac_fixed<W, I, false> &exp
  ){
    ac_fixed<W, I+1, 1> x_eff;
                        x_eff = x + (x >> 1) 
                                  - (x >> 4);

    ac_int<I-1, 0>      _int; 
    ac_fixed<W-I, 1, 1> _frac;
    split<W, I>(_int, _frac, x_eff);

    ac_fixed<W-I+1, 1, 0> _pow2;
    ac_fixed<W-I+1, 1, 0> a;
    ac_fixed<W-I+1, 1, 0> b;
    pow2_pwl<W, I>(a, b, _frac);
    
    _pow2 = a * _frac + b;
    
    exp = _pow2 >> _int;
  }

  template<int W=16, int I=5>
  void
  log2exp(
    const ac_fixed<W, I, true> &x,
          ac_int<I,false>      &out
  ){
    using dt = ac_fixed<W, I, true>;
    ac_fixed<W, I, true> x_eff = ~x+__AC_QNT(dt);
    x_eff = x_eff + (x_eff >> 1) - (x_eff >> 4);
    x_eff = x_eff[W-I-1] ? ac_fixed<W, I, true>(x_eff + 1) : x_eff;
    
    out = x_eff.template slc<I>(W-I);
  }

  template<int I>
  inline ac_int<CLOG2(I), false>
  lod(const ac_int<I, false> &x) {
    ac_int<CLOG2(I), false> lzc = reduction::lzcount<I>(x);
    return ac_int<CLOG2(I), false>(I - lzc - 1);
  }

  template<int W=16, int I=5>
  void
  aldivision(
    const ac_int<I,false>       &x,
          ac_fixed<W, I, false> &sum,
          ac_fixed<W, I, false> &out
  ){
    using dt = ac_fixed<W, I, false>;
    
    ac_int<I, false> lod_sum = lod(sum.template slc<I>(W-I));
    
    ac_int<I, false> diff = x + lod_sum;
    ac_fixed<W, I, false> shifted_sum = sum >> lod_sum;
    ac_int<1, false> sel_bit = shifted_sum.template slc<1>(W-I-1);
    dt seld_value = sel_bit ? dt(0.568f) : dt(0.818f);
    
    ac_int<W, false> a = seld_value.template slc<W>(0) >> (diff);

    out.set_slc(0, a);
  }

  // PWL approx of 2^x in (-1, 1)
  template< int W, int I >
  void pow2_pwl_sym(
          ac_fixed<W-I+1, 1, 0> &a,
          ac_fixed<W-I+1, 1, 0> &b,
    const ac_fixed<W-I+1, 1, 1> &x
  ){  
    /**
     * @note If x is zero, just output 1 no need to pwl.
    */

    const ac_fixed<W-I+1, 1, 0>
    alpha_lut[8] = {
      0.376789767154104,
      0.449069481221994,
      0.533857921384583,
      0.634687091204338,
      0.755034324221494,
      0.897990622883645,
      1.06640230597548,
      1.27408287473722
    };

    const ac_fixed<W-I+1, 1, 0> 
    beta_lut[8] = {
      0.87566839377251,
      0.929878179323428,
      0.972272399404722,
      0.997479691859661,
      0.997479691859661,
      0.961740617194123,
      0.877534775648206,
      0.721774349076904
    };

    /**
      * @param index is the 3 MSBs of the fraction.
      * Here x is <1>.<W-I> --> W-I+1 bits in total
    */
    ac_int<3, false> index;
                     index.set_slc(0, x.template slc<3>(W-I+1 - 4));
                     index[2] = ~x[W-I];
   
    ac_int<1, false> is_zero = (x == 0);

    a = is_zero ? ac_fixed<W-I+1, 1, 0>(0.0f) : alpha_lut[index];
    b = is_zero ? ac_fixed<W-I+1, 1, 0>(1.0f) : beta_lut[index] ;
  }
}

#endif