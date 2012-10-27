#!/usr/bin/python

import numpy
from numpy import *
from gnuradio import gr, window, digital
from gnuradio.extras import block_gateway
import logging


class gen_sss_fd:
  """
  Calculate the SSS in the frequency domain.
  """
  def __init__(self, N_id_1, N_id_2, N_re):
    logger = logging.getLogger('gen_sss_fd')
    
    self.N_id_1 = N_id_1
    self.N_id_2 = N_id_2
    (sss0, sss10) = self.generate(N_id_1, N_id_2)
    logger.debug(sss0)
    logger.debug(sss10)
    
    k = int((N_re)/2 - 31)
    self.sss0 = numpy.zeros(N_re,numpy.complex)
    self.sss10 = numpy.zeros(N_re,numpy.complex)
    self.sss0[k:k+62] = sss0
    self.sss10[k:k+62] = sss10


  def get_sss(self, slot0=True):
    if slot0:
      return self.sss0
    else:
      return self.sss10
    
    
  def get_sss_conj_rev(self, slot0=True):
    if slot0:
      return numpy.conjugate(self.sss0[::-1])
    else:
      return numpy.conjugate(self.sss10[::-1])

  
  def generate(self, N_id_1, N_id_2): 
    # Generate m0 and m1
    q_prime = int(floor(N_id_1/30))
    q       = int(floor((N_id_1 + (q_prime*(q_prime+1)/2))/30))
    m_prime = N_id_1 + (q*(q+1)/2)
    m0      = int(mod(m_prime, 31))
    m1      = int(mod((m0 + floor(m_prime/31) + 1), 31))

    # Generate s_tilda
    x_s_tilda = array([0,0,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1])
    s_tilda = 1-2*x_s_tilda

    # Generate c_tilda
    x_c_tilda = array([0,0,0,0,1,0,1,0,1,1,1,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0,1])
    c_tilda = 1-2*x_c_tilda

    # Generate z_tilda
    x_z_tilda = array([0,0,0,0,1,1,1,0,0,1,1,0,1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1])
    z_tilda = 1-2*x_z_tilda

    # Generate s0_m0 and s1_m1
    s0_m0 = s_tilda[mod(range(m0,31+m0),31)]
    s1_m1 = s_tilda[mod(range(m1,31+m1),31)]

    # Generate c0 and c1
    c0 = c_tilda[mod(range(N_id_2,31+N_id_2),31)]
    c1 = c_tilda[mod(range(N_id_2+3,31+N_id_2+3),31)]

    # Generate z1_m0 and z1_m1
    z1_m0 = z_tilda[mod(range(0,31)+mod(m0,8),31)]
    z1_m1 = z_tilda[mod(range(0,31)+mod(m1,8),31)]

    # Generate SSS
    sss_d_u_0 = numpy.zeros(62,numpy.complex)
    sss_d_u_0[0:62:2] = s0_m0 * c0
    sss_d_u_0[1:62:2] = s1_m1 * c1 * z1_m0 

    sss_d_u_10 = numpy.zeros(62,numpy.complex)
    sss_d_u_10[0:62:2] = s1_m1 * c0
    sss_d_u_10[1:62:2] = s0_m0 * c1 * z1_m1 

    return (sss_d_u_0,sss_d_u_10)
  
  

class gen_sss_td:
  def __init__(self, N_id_1, N_id_2, slot0=True, N_re=128, N_cp_ts=144, freq_corr=0):
    top = gr.top_block("foo");
    
    source = gr.vector_source_c(range(0,N_re), False, N_re)
    source.set_data(gen_sss_fd(N_id_1, N_id_2, N_re).get_sss(slot0)); 
    fft = gr.fft_vcc(N_re, False, window.blackmanharris(1024), True)
    cp = digital.ofdm_cyclic_prefixer(N_re, N_re+N_cp_ts*N_re/2048)
    if freq_corr != 0:
      freq_corr = gr.freq_xlating_fir_filter_ccf(1, 1, freq_corr, 15000*N_re)
    sink = gr.vector_sink_c(1)
    
    top.connect(source, fft, cp)
    if freq_corr != 0:
      top.connect(cp, freq_corr, sink)
    else:
      top.connect(cp, sink)
    top.run()
    self.data = sink.data()
    
    
  def get_data(self):
    return self.data
  
  def get_data_conj_rev(self):
    return numpy.conjugate(self.data[::-1]) 
