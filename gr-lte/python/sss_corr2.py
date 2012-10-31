#!/usr/bin/python

from numpy import zeros, argmax, ones
from gnuradio import gr, filter
from gnuradio.extras import block_gateway
from pss_source import *
from sss_source import *
from sss_equ2 import *  
from sss_derot import *
from sss_ml_fd2 import *
from symbol_source import *
import logging

class sss_corr2():
  """
  SSS correlator block
  """
  def __init__(self, decim=16, fft_size=2048/16, N_re=62, avg_frames=8, dump=None, N_id_1s = range(0,168), slot_0_10s = range(0,2)):
    self.logger = logging.getLogger('sss_corr2')
    
    # store parameters
    self.decim = decim
    self.fft_size = fft_size
    self.N_re = N_re
    self.avg_frames = avg_frames
    self.dump = dump
    self.N_id_1s = N_id_1s
    self.slot_0_10s = slot_0_10s

    # calculate statics    
    self.symbol_mask = numpy.zeros(20*7)
    self.symbol_mask[5:7] = 1
    self.symbol_mask[75:77] = 1

    # generate PSS sequences 
    self.pss_fd_vec = []
    for N_id_2 in range(0,3):
      self.pss_fd_vec.append(gen_pss_fd(N_id_2, self.N_re, False).get_data())
      
    # generate SSS sequences
    self.sss_fd_vec = []
    for N_id_2 in range(0,3):
      self.sss_fd_vec.append([])
      for N_id_1 in range(0,168):
        self.sss_fd_vec[N_id_2].append([])
        for slot_0_10 in range(0,2):
          self.sss_fd_vec[N_id_2][N_id_1].append(gen_sss_fd(N_id_1,N_id_2, N_re).get_sss_conj(slot_0_10!=0))
    self.pss_ref_src = gr.vector_source_c(zeros(0), True, self.N_re)
    self.sss_ref_src = gr.vector_source_c(zeros(0), True, self.N_re)

    # SSS equalization flow graph
    self.equ_source = symbol_source(decim=self.decim, vlen=self.fft_size)
    fft = gr.fft_vcc(self.fft_size, True, (window.blackmanharris(1024)), True, 1)
    vector_to_stream_0 = gr.vector_to_stream(gr.sizeof_gr_complex*1, self.fft_size)
    keep_m_in_n_0 = gr.keep_m_in_n(gr.sizeof_gr_complex, self.N_re/2, self.fft_size, (self.fft_size-self.N_re)/2-1)
    keep_m_in_n_1 = gr.keep_m_in_n(gr.sizeof_gr_complex, self.N_re/2, self.fft_size, (self.fft_size)/2)
    stream_mux_0 = gr.stream_mux(gr.sizeof_gr_complex, 2*[self.N_re/2]) 
    stream_to_vector_0_0 = gr.stream_to_vector(gr.sizeof_gr_complex*1, self.N_re)
    deinterleave_0 = gr.deinterleave(gr.sizeof_gr_complex*self.N_re)
    self.equ = sss_equ2()
    self.equ_sink = gr.vector_sink_c(self.N_re)

    self.equ_top = gr.top_block("sss equ graph")
    self.equ_top.connect(self.equ_source, fft, vector_to_stream_0)
    self.equ_top.connect(vector_to_stream_0, keep_m_in_n_0, (stream_mux_0,0))
    self.equ_top.connect(vector_to_stream_0, keep_m_in_n_1, (stream_mux_0,1))
    self.equ_top.connect(stream_mux_0, stream_to_vector_0_0, deinterleave_0)
    self.equ_top.connect((deinterleave_0,0), (self.equ,0))
    self.equ_top.connect((deinterleave_0,1), (self.equ,1))
    self.equ_top.connect(self.pss_ref_src, (self.equ,2))
    self.equ_top.connect((self.equ,0), self.equ_sink)
    self.equ_top.connect((self.equ,1), gr.null_sink(self.N_re*gr.sizeof_gr_complex))

    if self.dump != None:
      self.equ_top.connect(self.equ_source, gr.file_sink(gr.sizeof_gr_complex*self.fft_size, self.dump + "_sss_td_in.cfile"))
      self.equ_top.connect((deinterleave_0,0), gr.file_sink(gr.sizeof_gr_complex*self.N_re, self.dump + "_sss_fd_in.cfile"))
      self.equ_top.connect((deinterleave_0,1), gr.file_sink(gr.sizeof_gr_complex*self.N_re, self.dump + "_pss_fd_in.cfile"))
      self.equ_top.connect(self.equ, gr.file_sink(gr.sizeof_gr_complex*self.N_re, self.dump + "_sss_fd_equ.cfile"))
    
    # SSS maximum likelihood estimation 
    self.ml_src = gr.vector_source_c(zeros(0), False, self.N_re)
    self.ml_derot = sss_derot()
    self.ml_sss = sss_ml_fd2(avg_frames=self.avg_frames)
    self.ml_sink = gr.vector_sink_f()

    self.ml_top = gr.top_block("sss ml graph")
    self.ml_top.connect(self.ml_src, self.ml_derot, self.ml_sss, self.ml_sink)
    self.ml_top.connect(self.sss_ref_src, (self.ml_derot,1))
    self.ml_top.connect(self.sss_ref_src, (self.ml_sss,1))

    if self.dump != None:
      self.ml_top.connect(self.ml_derot, gr.file_sink(gr.sizeof_gr_complex*self.N_re, self.dump + "_sss_fd_derot.cfile"))
      self.ml_top.connect(self.ml_sss, gr.file_sink(gr.sizeof_float, self.dump + "_sss_corr.cfile"))


  def run_equ(self, buffer, frame_time, N_id_2):
    # setup
    self.equ_source.set_data(buffer, self.symbol_mask, frame_time)
    self.pss_ref_src.set_data(self.pss_fd_vec[N_id_2])
    self.equ_sink.reset()
    # run equalization
    self.equ_top.run()


  def run_ml(self, N_id_1, N_id_2, slot_0_10):
    # prepare
    self.ml_src.set_data(self.equ_sink.data())
    self.sss_ref_src.set_data(concatenate((self.sss_fd_vec[N_id_2][N_id_1][slot_0_10==1],self.sss_fd_vec[N_id_2][N_id_1][slot_0_10!=1]),1))
    self.ml_sink.reset()
    # run
    self.ml_top.run()
    # collect result
    ml_data = self.ml_sink.data()
    return mean(ml_data)


  def correlate(self, buffer, frame_time, N_id_2):
    self.run_equ(buffer, frame_time, N_id_2)
    ml_sss = -1e5*numpy.ones(168*2)
    for N_id_1 in self.N_id_1s:
      for slot_0_10 in self.slot_0_10s:
        ml_sss[2*N_id_1+slot_0_10] = self.run_ml(N_id_1, N_id_2, slot_0_10)
    self.logger.debug("SSS ML res\n" + str(ml_sss))
    idx = argmax(ml_sss)
    if idx % 2 > 0:
      frame_time += 30720*5
    return [(int(idx/2), frame_time, ml_sss[idx])]





