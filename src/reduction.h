#ifndef __REDUCTION_H__
#define __REDUCTION_H__

namespace reduction{
  template<int N>
  struct max_s {
  template<typename T>
    static T max(T *a) {
      T m0 = max_s<N/2>::max(a);
      T m1 = max_s<N-N/2>::max(a+N/2);
      
      return m0 > m1 ? m0 : m1;
    }
  };

  template<> 
  struct max_s<1> {
    template<typename T>
    static T max(T *a) {
      return a[0];
    }
  };

  template<int N, typename T>
  T max(T *a) {
    return max_s<N>::max(a);
  };

  template <int W>
  bool lzc_reduce (ac_int<W,false> a) {
    ac_int<W/2,false> odd, even;

    // split odd even bits of a
    #pragma unroll yes
    SPLIT:for (int i=0; i < W/2; i++) {
      even[i] = a[2*i];
      odd[i] = a[2*i+1];
    }

    // prefix AND from MSB-to-LSB of inverted even bits
    even = even.bit_complement();
    ac_int<W/2,false> even_prefix;
    ac_int<1,false> t=true;
    
    #pragma unroll yes
    PREFIXAND:for (int i=W/2-1; i >=0; i--) {
      if (i == W/2-1) t = even[i];
      else t = t & even[i];
      even_prefix[i] = t;
    }

    // fix alignment of prefixed and terms
    even_prefix = even_prefix >> 1;
    even_prefix[W/2-1] = 1;

    // prepare terms for each bit position
    ac_int<W/2,false> tmp = even_prefix & odd;

    // return the wide OR of those terms
    return tmp.or_reduce();
  }

  // Version 1: try to do it
  template<int N>
  struct lzc_s {
    template<typename T>
    static void lzc(ac_int<N,false> a, T &out) {
      ac_int<N/2, false> a0;
      out[ac::log2_ceil<N>::val] = lzc_reduce<N>(a);

      #pragma unroll yes
      for (int i = 0; i < N / 2; i++) {
        a0[i] = a[2*i] | a[2*i + 1];
      }
      lzc_s<N/2>::lzc(a0,out);
    }
  };

  template<>
  struct lzc_s<1> {
    template<typename T>
    static void lzc(ac_int<1,false> a,  T &out) {
      out[0] = a[0];
    }
  };


  template<int N=8>
  ac_int<ac::log2_ceil<N>::val+1,false> lzcount(ac_int<N,false> x) {
    ac_int<ac::log2_ceil<N>::val+1,false> b,res;
    lzc_s<N>::lzc(x,b);

    // reverse bits
    #pragma unroll yes
    for (int i=0; i< ac::log2_ceil<N>::val+1; i++) {
      res[i] = b[ac::log2_ceil<N>::val - i]; 
    }

    // complement them and send them out
    return (res.or_reduce() == 0) ? (ac_int<ac::log2_ceil<N>::val+1,false>)0 :  res.bit_complement();
  }
}

#endif