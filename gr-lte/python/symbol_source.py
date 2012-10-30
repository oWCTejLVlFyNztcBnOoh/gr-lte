#!/usr/bin/python

import numpy
from numpy import zeros, mean
from gnuradio import gr, window, digital
from gnuradio.extras import block_gateway
import logging

class symbol_source(gr.hier_block2):
    """
    Cut symbols from sample stream and put them to output 
    @todo: implement extended CP. currently, normal CP is assumed.
    @todo: implement del_cp = False 
    """
    
    def __init__(self, data=None, decim=16, symbol_mask=None, frame_start=-1, del_cp=True, symbol_start=-1, repeat=False, vlen=1):
      """
      @param decim: sample rate decimator, sample rate = 30720e3/decim
      @param symbol_mask: array[0..10*14] in {0,1} of symbols to deliver per frame
      @param frame_start: start of the frame relative to data in Ts
      @param del_cp: whether to remove the CP    
      @param symbol_start: start position of the FFT window within (CP+symbol)
      @param repeat: whether to repeat after data is finished 
      """
      gr.hier_block2.__init__(
          self, "symbol source",
          gr.io_signature(0, 0, 0),
          gr.io_signature(1, 1, vlen*gr.sizeof_gr_complex),
      )
      self.logger = logging.getLogger('symbol_source')
      
      self.vlen = vlen
      self.repeat = repeat
      self.decim = decim

      samp_cp_n = int(144/self.decim)
      if symbol_start < 0:
        self.symbol_start = samp_cp_n
      else:
        self.symbol_start = int(symbol_start/self.decim)
      
      self.source = gr.vector_source_c(zeros(0), self.repeat, self.vlen)
      self.connect(self.source, self)
      if None != data:
        self.set_data(data, symbol_mask, frame_start)

    
    def set_data(self, data, symbol_mask, frame_start):
      samp_per_symb = int(2048/self.decim)
      samp_cp_n = int(144/self.decim)
      samp_cp_e = int(160/self.decim)

      pos_i = int(frame_start/self.decim)
      symbols = zeros((len(data)-pos_i)*mean(symbol_mask), numpy.complex)
      pos_o = 0
      n_symb_per_slot = 0
      n_slot = 0
      while pos_i+samp_per_symb+samp_cp_e < len(data):
        samp_cp = samp_cp_n
        if 0 == n_symb_per_slot:
          samp_cp = samp_cp_e
        if symbol_mask[n_slot*7+n_symb_per_slot] > 0:
          osta = pos_o
          oend = pos_o+samp_per_symb
          ista = pos_i+self.symbol_start
          iend = pos_i+self.symbol_start+samp_per_symb
          self.logger.debug("Cut slot/symbol {:2d}/{:1d} from in[{}:{}] to out[{}:{}]".format(n_slot,n_symb_per_slot,ista,iend,osta,oend))
          symbols[pos_o:pos_o+samp_per_symb] = data[pos_i+self.symbol_start:pos_i+self.symbol_start+samp_per_symb]
          pos_o += samp_per_symb
        pos_i += samp_per_symb + samp_cp
        n_symb_per_slot += 1
        if n_symb_per_slot >= 7:
          n_symb_per_slot = 0
          n_slot += 1
          if n_slot >= 20:
            n_slot = 0
      symbols = symbols[0:pos_o]
      
      self.source.set_data(symbols)
      self.source.rewind()

    