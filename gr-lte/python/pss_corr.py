#!/usr/bin/python

import numpy
from gnuradio import gr, filter
from gnuradio.extras import block_gateway
from pss_source import *
from sss_source import *

class pss_corr(gr.hier_block2):
  """
  PSS correlator block
  """
  def __init__(self, N_id_2, decim=16, avg_halfframes=2*8, freq_corr=0):
    gr.hier_block2.__init__(
        self, "PSS correlator",
        gr.io_signature(1, 1, gr.sizeof_gr_complex),
        gr.io_signature(1, 1, gr.sizeof_int),
    )
    vec_half_frame = 30720*5/decim
    
    taps = gen_pss_td(N_id_2, N_re=2048/decim, freq_corr=freq_corr).get_data_conj_rev()
    self.corr = filter.fir_filter_ccc(1, taps)
    self.mag = gr.complex_to_mag_squared()
    self.vec = gr.stream_to_vector(gr.sizeof_float*1, vec_half_frame)
    self.deint = gr.deinterleave(gr.sizeof_float*vec_half_frame)
    self.add = gr.add_vff(vec_half_frame)
    self.argmax = gr.argmax_fs(vec_half_frame)
    self.null = gr.null_sink(gr.sizeof_short*1)
    self.to_float = gr.short_to_float(1, 1)
    self.to_int = gr.float_to_int(1, decim)
    #self.framestart = gr.add_const_ii(-160-144*5-2048*6+30720*5)
    
    self.connect(self, self.corr, self.mag, self.vec)
    self.connect((self.argmax,1), self.null)
    #self.connect(self.argmax, self.to_float, self.to_int, self.framestart, self)
    self.connect(self.argmax, self.to_float, self.to_int, self)
    
    if avg_halfframes == 1:
      self.connect(self.vec, self.argmax)
    else:
      self.connect(self.vec, self.deint)
      self.connect(self.add, self.argmax)
      for i in range(0, avg_halfframes):
        self.connect((self.deint, i), (self.add, i))







