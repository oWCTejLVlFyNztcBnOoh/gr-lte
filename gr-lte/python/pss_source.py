#!/usr/bin/python

import numpy
from numpy import cos, sin
from gnuradio import gr, window, digital
from gnuradio.extras import block_gateway

class pss_source_fd(gr.hier_block2):
    """
    Creates a PSS symbol in the frequency domain
    Document Reference: 3GPP TS 36.211 v10.1.0 section 6.11.1.1
    """
    def __init__(self, N_id_2, N_re=128, repeat=False):
      gr.hier_block2.__init__(
          self, "PSS source frequency-domain",
          gr.io_signature(0, 0, 0),
          gr.io_signature(1, 1, gr.sizeof_gr_complex*N_re),
      )
      self.N_re = N_re

      pss = self.generate_pss(N_id_2)
      pss_fd = numpy.zeros(self.N_re,numpy.complex)
      pss_fd[1:32] = pss[31:62]
      pss_fd[self.N_re-31:] = pss[0:31]

      self.source = gr.vector_source_c(range(0,self.N_re), repeat, self.N_re)
      self.source.set_data(pss_fd);
      self.connect(self.source, self)
      
    
    def generate_pss(self,N_id_2):
        if N_id_2 == 0:
            root_idx = 25
        elif N_id_2 == 1:
            root_idx = 29
        else:
            root_idx = 34
    
        pss = numpy.zeros(62,numpy.complex)
        for i in range(0,31):
            pss[i] = numpy.complex(cos(-numpy.pi*root_idx*i*(i+1)/63), sin(-numpy.pi*root_idx*i*(i+1)/63))
        for i in range(31,62):
            pss[i] = numpy.complex(cos(-numpy.pi*root_idx*(i+1)*(i+2)/63), sin(-numpy.pi*root_idx*(i+1)*(i+2)/63))
            
        return (pss)


class pss_source_td(gr.hier_block2):
  """
  Creates a PSS symbol in the time domain, adds cyclic prefix + frequency correction (optional)
  Document Reference: 3GPP TS 36.211 v10.1.0 section 6.11.1.1
  """
  def __init__(self, N_id_2, N_re=128, N_cp_ts=144, freq_corr=0, repeat=False):
    gr.hier_block2.__init__(
        self, "PSS source time-domain",
        gr.io_signature(0, 0, 0),
        gr.io_signature(1, 1, gr.sizeof_gr_complex),
    )
    
    self.pss_fd = pss_source_fd(N_id_2, N_re, repeat);
    self.fft = gr.fft_vcc(N_re, False, window.blackmanharris(1024), False)
    self.cp = digital.ofdm_cyclic_prefixer(N_re, N_re+N_cp_ts*N_re/2048)
    if freq_corr != 0:
      self.freq_corr = gr.freq_xlating_fir_filter_ccf(1, 1, freq_corr, 15000*N_re)
    
    self.connect(self.pss_fd, self.fft, self.cp)
    if freq_corr != 0:
      self.connect(self.cp, self.freq_corr, self)
    else:
      self.connect(self.cp, self)


def gen_pss_td(N_id_2, N_re=128, N_cp_ts=144, freq_corr=0):
  top = gr.top_block("foo");
  pss_src_td = pss_source_td(N_id_2, N_re, N_cp_ts, freq_corr)
  sink = gr.vector_sink_c(1)
  top.connect(pss_src_td, sink)
  top.run()
  return sink.data()
  

class gen_pss_td:
  def __init__(self,N_id_2, N_re=128, N_cp_ts=144, freq_corr=0):
    top = gr.top_block("foo");
    pss_src_td = pss_source_td(N_id_2, N_re, N_cp_ts, freq_corr)
    sink = gr.vector_sink_c(1)
    top.connect(pss_src_td, sink)
    top.run()
    self.data = sink.data()
    
  def get_data(self):
    return self.data
  
  def get_data_conj_rev(self):
    return numpy.conjugate(self.data[::-1]) 







