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
  def __init__(self, N_id_2, decim=16, avg_halfframes=2*8, freq_corr=0, dump=None):
    gr.hier_block2.__init__(
        self, "PSS correlator",
        gr.io_signature(1, 1, gr.sizeof_gr_complex),
        gr.io_signature(1, 1, gr.sizeof_float),
    )
    vec_half_frame = 30720*5/decim
    
    self.taps = []
    for i in range(0,3):
      self.taps.append(gen_pss_td(i, N_re=2048/decim, freq_corr=freq_corr).get_data_conj_rev())
    self.corr = filter.fir_filter_ccc(1, self.taps[N_id_2])
    self.mag = gr.complex_to_mag_squared()
    self.vec = gr.stream_to_vector(gr.sizeof_float*1, vec_half_frame)
    self.deint = gr.deinterleave(gr.sizeof_float*vec_half_frame)
    self.add = gr.add_vff(vec_half_frame)
    self.argmax = gr.argmax_fs(vec_half_frame)
    self.null = gr.null_sink(gr.sizeof_short*1)
    self.max = gr.max_ff(vec_half_frame)
    self.to_float = gr.short_to_float(1, 1./decim)
    self.interleave = gr.interleave(gr.sizeof_float)
    #self.framestart = gr.add_const_ii(-160-144*5-2048*6+30720*5)
    
    self.connect(self, self.corr, self.mag, self.vec)
    self.connect((self.argmax,1), self.null)
    #self.connect(self.argmax, self.to_float, self.to_int, self.framestart, self)
    self.connect(self.argmax, self.to_float, self.interleave, self)
    self.connect(self.max, (self.interleave,1))
    
    if avg_halfframes == 1:
      self.connect(self.vec, self.argmax)
      self.connect(self.vec, self.max)
    else:
      self.connect(self.vec, self.deint)
      self.connect(self.add, self.argmax)
      self.connect(self.add, self.max)
      for i in range(0, avg_halfframes):
        self.connect((self.deint, i), (self.add, i))

    if dump != None:
      self.connect(self.mag, gr.file_sink(gr.sizeof_float, dump + "_pss_corr_f.cfile"))
      self.connect(self.add, gr.file_sink(gr.sizeof_float*vec_half_frame, dump + "_pss_corr_add_f.cfile"))


  def set_N_id_2(self, N_id_2):
    self.corr.set_taps(self.taps[N_id_2])




