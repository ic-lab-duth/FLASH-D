#ifndef __FLATTEN_H__
#define __FLATTEN_H__

#include "defines.h"
#include "fp_arithm.h"
#include "ac_channel.h"

#ifndef __SYNTHESIS__
#include <math.h>
#endif

#include "mc_scverify.h"

template <int d=64>
void fa2_dotProd(
  ac_channel<vec_t>  &Q_i,
  ac_channel<vec_t>  &K_i,
  ac_channel<bit_t > &in2dot_i,
  ac_channel<fptype> &s_o,
  ac_channel<bit_t > &dot2ol_o
){
  vec_t q_i, k_i;
  bit_t is_last;

  #ifndef __SYNTHESIS__
  if(
    Q_i.available(1) && 
    K_i.available(1) &&
    in2dot_i.available(1)
  )
  #endif
  {
    Q_i.read(q_i);
    K_i.read(k_i);
    in2dot_i.read(is_last);

    fptype x = m_dotProd<fptype, ffptype, d>(q_i.data, k_i.data);
    
    s_o.write(x);
    dot2ol_o.write(is_last);
  }
}

template <int d=64>
void CCS_BLOCK(fa2_online_acc)(
  ac_channel<fptype> &s_i,
  ac_channel<bit_t > &dot2ol_i,
  ac_channel<fptype> &corr_o,
  ac_channel<fptype> &e_xmax_o,
  ac_channel<fptype> &l_o,
  ac_channel<bit_t > &ol2div_o
){
  static fptype max_old = fptype(-1000.0f);
  static fptype sum_old = fptype(0.0f);

  bit_t is_last;
  fptype x;
  #ifndef __SYNTHESIS__
  if(
    s_i.available(1) &&
    dot2ol_i.available(1)
  )
  #endif
  {
    s_i.read(x);
    dot2ol_i.read(is_last);

    fptype max_x = (x > max_old) ? x : max_old;
    fptype norm_x   = x - max_x;
    fptype diff_max = max_old - max_x;
    
    fptype exp_x   = m_ac_exp_pwl(norm_x);
    fptype exp_max = m_ac_exp_pwl(diff_max);
    
    fptype sum_old_x_exp_max = sum_old * exp_max;
    fptype sum               = sum_old_x_exp_max + exp_x;

    max_old = max_x;
    sum_old = sum;

    corr_o  .write(exp_max);
    e_xmax_o.write(exp_x);
    l_o     .write(sum);
    ol2div_o.write(is_last);
  }
}

template <int d=64>
void CCS_BLOCK(fa2_output_acc)(
  ac_channel<fptype> &corr_i,  
  ac_channel<fptype> &e_xmax_i,
  ac_channel<vec_t>  &V_i,
  ac_channel<bit_t>  &in2oacc_i,
  ac_channel<vec_t>  &O_o,
  ac_channel<bit_t>  &oacc2div_o
){
  static vec_t o_old;

  bit_t is_last;
  fptype corr, e_xmax;
  vec_t v_i;
  #ifndef __SYNTHESIS__
  if(
    V_i      .available(1) &&
    corr_i   .available(1) &&
    e_xmax_i .available(1) &&
    in2oacc_i.available(1)
  )
  #endif
  {
    in2oacc_i.read(is_last);
    corr_i   .read(corr);
    e_xmax_i .read(e_xmax);
    V_i      .read(v_i);
    
    #pragma hls_unroll yes
    O_ACC: for(int i = 0; i < d; i++) {
      o_old.data[i] = o_old.data[i] * corr + e_xmax * v_i.data[i];
    }
    
    O_o.write(o_old);
    oacc2div_o.write(is_last);
  }
}

template <int d=64>
void CCS_BLOCK(fa2_div)(
  ac_channel<vec_t>  &o_i,
  ac_channel<fptype> &l_i,
  ac_channel<bit_t>  &oacc2div_i,
  ac_channel<bit_t>  &ol2div_i,
  ac_channel<vec_t>  &O_o
){
  bit_t is_last_ol, is_last_oa;
  vec_t  o;
  fptype l;
  vec_t  o_out;

  #ifndef __SYNTHESIS__
  if(
    o_i.available(1)        && 
    l_i.available(1)        &&
    oacc2div_i.available(1) &&
    ol2div_i.available(1)
  )
  #endif
  {
    o_i       .read(o);
    l_i       .read(l);
    oacc2div_i.read(is_last_oa);
    ol2div_i  .read(is_last_ol);

    #pragma hls_unroll yes
    DIV: for(int j = 0; j < d; j++) {
      o_out.data[j] = o.data[j] * m_ac_recip_pwl(l);
    }

    if(is_last_oa && is_last_ol) {
      O_o.write(o_out);
    }
  }
}

#pragma hls_design top
template <int Tc, int Bc, int d=64>
void CCS_BLOCK(flashAtten2)(
  fptype Q_i[d],
  fptype K_i[Tc*Bc][d],
  fptype V_i[Tc*Bc][d],

  fptype O_o[d]
){
  static ac_channel<vec_t> Q_ch;
  static ac_channel<vec_t> K_ch;
  static ac_channel<fptype> s_ch;
  
  static ac_channel<fptype> corr_ch;
  static ac_channel<fptype> e_xmax_ch;
  static ac_channel<fptype> l_ch;
  
  static ac_channel<vec_t> V_ch;
  static ac_channel<vec_t> o_ch;
  
  static ac_channel<vec_t> O_ch;

  // `is_last` control signal
  static ac_channel<bit_t> in2dot_ch;
  static ac_channel<bit_t> in2oa_ch;
  static ac_channel<bit_t> dot2ol_ch;
  static ac_channel<bit_t> ol2div_ch;
  static ac_channel<bit_t> oacc2div_ch;

  vec_t q;
  #pragma hls_unroll yes
  WRITE_Q: for(int j = 0; j < d; j++) {
    q.data[j] = Q_i[j];
  }
  
  vec_t k, v;
  #pragma hls_pipeline_init_interval 1
  WRITE_KV: for(int i = 0; i < Tc*Bc; i++) {
    #pragma hls_unroll yes
    WRITE_KV_i: for(int j = 0; j < d; j++) {
      k.data[j] = K_i[i][j];
      v.data[j] = V_i[i][j];
    }
    K_ch.write(k);
    V_ch.write(v);
    Q_ch.write(q);

    in2dot_ch.write(bit_t(i==(Tc*Bc-1)));
    in2oa_ch.write(bit_t(i==(Tc*Bc-1)));
  }

  #ifndef __SYNTHESIS__
  while(
    K_ch.available(1) &&
    V_ch.available(1)
  )
  #endif
  {
    fa2_dotProd   <d>(Q_ch   , K_ch     , in2dot_ch  , s_ch      , dot2ol_ch         );
    fa2_online_acc<d>(s_ch   , dot2ol_ch, corr_ch    , e_xmax_ch , l_ch , ol2div_ch  );
    fa2_output_acc<d>(corr_ch, e_xmax_ch, V_ch       , in2oa_ch , o_ch , oacc2div_ch );
    fa2_div       <d>(o_ch   , l_ch     , oacc2div_ch, ol2div_ch , O_ch              );
  }

  vec_t o;
  #ifndef __SYNTHESIS
  while(O_ch.available(1))
  #endif
  {
    O_ch.read(o);
  }

  #pragma hls_unroll yes
  WRITE_O: for(int i = 0; i < d; i++) {
    O_o[i] = o.data[i];
  }
}

template <int d=64>
void CCS_BLOCK(w_dotProd)(
  ac_channel<vec_t>  &Q_i,
  ac_channel<vec_t>  &K_i,
  ac_channel<fptype> &s_o
){
  vec_t q_i, k_i;

  #ifndef __SYNTHESIS__
  if(
    Q_i.available(1) && 
    K_i.available(1)
  )
  #endif
  {
    Q_i.read(q_i);
    K_i.read(k_i);
    
    fptype x = m_dotProd<fptype, ffptype, d>(q_i.data, k_i.data);
    
    s_o.write(x);
  }
}


template <int d=64>
void CCS_BLOCK(w_online_acc)(
  ac_channel<fptype>           &s_i,
  ac_channel<fptype>           &w_o,
  ac_channel<ac_int<3,false> > &onehot_o
){
  static fptype x_old = fptype(0.0f);
  static fptype w_old = fptype(1.0f);

  ac_int<2, false> one_hot;
  fptype buff_w;
  fptype x;
  #ifndef __SYNTHESIS__
  if(
    s_i.available(1)
  )
  #endif
  {
    s_i.read(x);

    fptype dx = x - x_old;

    one_hot[0] = (dx <= fptype(-6.0f));
    one_hot[1] = (dx >= fptype(11.0f));

    fptype w  = one_hot[0] ? fptype(0.00001f) : 
                one_hot[1] ? fptype(0.99999f) : m_sigmoid_pwl(dx + m_ln_pwl(w_old));

    x_old = x;
    w_old = w; 

    w_o      .write(w);
    onehot_o .write(one_hot);
  }
}

template <int d=64>
void CCS_BLOCK(w_output_acc)(
  ac_channel<fptype>           &w_i       ,
  ac_channel<vec_t >           &V_i       ,
  ac_channel<ac_int<3,false> > &onehot_i  ,
  ac_channel<bit_t >           &in2oacc_i ,
  ac_channel<vec_t >           &O_o       
){
  static vec_t o_old;

  ac_int<3,false> one_hot;
  bit_t is_last;
  fptype w;
  vec_t v_i;
  #ifndef __SYNTHESIS__
  if(
    w_i.available(1) && 
    V_i.available(1) && 
    in2oacc_i.available(1) &&
    onehot_i.available(1)
  )
  #endif
  {
    V_i.read(v_i);
    onehot_i.read(one_hot);
    in2oacc_i.read(is_last);
    w_i.read(w);

    #pragma hls_unroll yes
    O_ACC: for (int i = 0; i < d; i++) {
      switch(one_hot) {
        case 1: {
          o_old.data[i] = o_old.data[i];
          break;
        }
        
        case 2: {
          o_old.data[i] = v_i.data[i];
          break;
        }
  
        default: {
          o_old.data[i] = w * v_i.data[i] + (fptype(1)-w)*o_old.data[i];
          break;
        }
      }
    }
    
    if (is_last) {
      O_o.write(o_old);
    }
  }
}

#pragma hls_design top
template <int Tc, int Bc, int d=64>
void CCS_BLOCK(weightAtten2)(
  fptype Q_i[d],
  fptype K_i[Tc*Bc][d],
  fptype V_i[Tc*Bc][d],

  fptype O_o[d]
){
  typedef ac_int<3,false> sel_t;

  static ac_channel<vec_t > Q_ch;
  static ac_channel<vec_t > K_ch;
  static ac_channel<fptype> s_ch;
  
  static ac_channel<fptype> w_ch;
  static ac_channel<sel_t > onehot_ch;

  static ac_channel<vec_t> V_ch;
  static ac_channel<vec_t> O_ch;

  // `is_last` control signal
  static ac_channel<bit_t> in2dot_ch;
  static ac_channel<bit_t> in2oa_ch;

  vec_t q;
  #pragma hls_unroll yes
  WRITE_Q: for(int j = 0; j < d; j++) {
    q.data[j] = Q_i[j];
  }
  
  vec_t k, v;
  WRITE_KV: for(int i = 0; i < Tc*Bc; i++) {
    #pragma hls_unroll yes
    WRITE_KV_i: for(int j = 0; j < d; j++) {
      k.data[j] = K_i[i][j];
      v.data[j] = V_i[i][j];
    }
    K_ch.write(k);
    V_ch.write(v);
    Q_ch.write(q);

    in2oa_ch.write(bit_t(i==(Tc*Bc-1)));
  }

  #ifndef __SYNTHESIS__
  while(
    K_ch.available(1) &&
    V_ch.available(1)
  )
  #endif
  {
    w_dotProd   <d>(Q_ch, K_ch , s_ch);
    w_online_acc<d>(s_ch, w_ch , onehot_ch);
    w_output_acc<d>(w_ch, V_ch , onehot_ch, in2oa_ch, O_ch);
  }

  vec_t o;
  #ifndef __SYNTHESIS__
  while(O_ch.available(1))
  #endif
  {
    O_ch.read(o);
  }

  #pragma hls_unroll yes
  WRITE_O: for(int i = 0; i < d; i++) {
    O_o[i] = o.data[i];
  }
}

#endif