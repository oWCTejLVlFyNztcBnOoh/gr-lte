#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: LTE DL synchronization signal
# Generated: Sat Oct 27 18:20:11 2012
##################################################

from gnuradio import digital
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import window
from gnuradio.eng_option import eng_option
from gnuradio.gr import firdes
from gnuradio.wxgui import scopesink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
from pss_corr import *  
from sss_corr import *   
import wx

class lte_dl_ss_source(grc_wxgui.top_block_gui):

	def __init__(self, N_id_1=134, N_id_2=0, decim=16, frames=1, gain=.08):
		grc_wxgui.top_block_gui.__init__(self, title="LTE DL synchronization signal")
		_icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
		self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

		##################################################
		# Parameters
		##################################################
		self.N_id_1 = N_id_1
		self.N_id_2 = N_id_2
		self.decim = decim
		self.frames = frames
		self.gain = gain

		##################################################
		# Variables
		##################################################
		self.samp_rate = samp_rate = 30720e3/decim
		self.fft_size = fft_size = 2048/decim
		self.N_re = N_re = 2048/decim

		##################################################
		# Blocks
		##################################################
		self.wxgui_scopesink2_0_2 = scopesink2.scope_sink_c(
			self.GetWin(),
			title="Scope Plot",
			sample_rate=samp_rate,
			v_scale=0,
			v_offset=0,
			t_scale=0,
			ac_couple=False,
			xy_mode=False,
			num_inputs=1,
			trig_mode=gr.gr_TRIG_MODE_AUTO,
			y_axis_label="Counts",
		)
		self.Add(self.wxgui_scopesink2_0_2.win)
		self.wxgui_scopesink2_0 = scopesink2.scope_sink_c(
			self.GetWin(),
			title="Scope Plot",
			sample_rate=samp_rate,
			v_scale=0,
			v_offset=0,
			t_scale=0,
			ac_couple=False,
			xy_mode=False,
			num_inputs=1,
			trig_mode=gr.gr_TRIG_MODE_AUTO,
			y_axis_label="Counts",
		)
		self.Add(self.wxgui_scopesink2_0.win)
		self.gr_vector_to_stream_1_0_0 = gr.vector_to_stream(gr.sizeof_gr_complex*1, fft_size+9)
		self.gr_vector_source_x_0_0_2_0_0 = gr.vector_source_c((gen_pss_td(N_id_2, N_re, 144).get_data()), True, fft_size+9)
		self.gr_vector_source_x_0_0_2_0 = gr.vector_source_c((gen_sss_td(N_id_1, N_id_2, True, N_re).get_data()), True, fft_size+9)
		self.gr_vector_source_x_0_0_0 = gr.vector_source_c((gen_pss_fd(N_id_2, fft_size, False).get_data()), True, fft_size)
		self.gr_vector_source_x_0 = gr.vector_source_c((gen_sss_fd( N_id_1, N_id_2, fft_size).get_sss(True)), True, fft_size)
		self.gr_throttle_0_0 = gr.throttle(gr.sizeof_gr_complex*1, samp_rate)
		self.gr_throttle_0 = gr.throttle(gr.sizeof_gr_complex*1, samp_rate)
		self.gr_null_sink_1 = gr.null_sink(gr.sizeof_gr_complex*1)
		self.gr_null_sink_0 = gr.null_sink(gr.sizeof_gr_complex*1)
		self.gr_multiply_const_vxx_0_0 = gr.multiply_const_vcc((1, ))
		self.gr_interleave_0_0 = gr.interleave(gr.sizeof_gr_complex*137)
		self.gr_interleave_0 = gr.interleave(gr.sizeof_gr_complex*fft_size)
		self.gr_fft_vxx_1 = gr.fft_vcc(fft_size, False, (window.blackmanharris(1024)), True, 1)
		self.digital_ofdm_cyclic_prefixer_0 = digital.ofdm_cyclic_prefixer(fft_size, fft_size+144/decim)

		##################################################
		# Connections
		##################################################
		self.connect((self.gr_interleave_0, 0), (self.gr_fft_vxx_1, 0))
		self.connect((self.gr_fft_vxx_1, 0), (self.digital_ofdm_cyclic_prefixer_0, 0))
		self.connect((self.gr_vector_source_x_0, 0), (self.gr_interleave_0, 1))
		self.connect((self.gr_vector_source_x_0_0_0, 0), (self.gr_interleave_0, 0))
		self.connect((self.gr_multiply_const_vxx_0_0, 0), (self.gr_null_sink_0, 0))
		self.connect((self.digital_ofdm_cyclic_prefixer_0, 0), (self.gr_null_sink_1, 0))
		self.connect((self.gr_vector_source_x_0_0_2_0, 0), (self.gr_interleave_0_0, 0))
		self.connect((self.gr_vector_source_x_0_0_2_0_0, 0), (self.gr_interleave_0_0, 1))
		self.connect((self.gr_interleave_0_0, 0), (self.gr_vector_to_stream_1_0_0, 0))
		self.connect((self.gr_vector_to_stream_1_0_0, 0), (self.gr_multiply_const_vxx_0_0, 0))
		self.connect((self.digital_ofdm_cyclic_prefixer_0, 0), (self.gr_throttle_0_0, 0))
		self.connect((self.gr_throttle_0_0, 0), (self.wxgui_scopesink2_0, 0))
		self.connect((self.gr_multiply_const_vxx_0_0, 0), (self.gr_throttle_0, 0))
		self.connect((self.gr_throttle_0, 0), (self.wxgui_scopesink2_0_2, 0))

	def get_N_id_1(self):
		return self.N_id_1

	def set_N_id_1(self, N_id_1):
		self.N_id_1 = N_id_1

	def get_N_id_2(self):
		return self.N_id_2

	def set_N_id_2(self, N_id_2):
		self.N_id_2 = N_id_2

	def get_decim(self):
		return self.decim

	def set_decim(self, decim):
		self.decim = decim
		self.set_N_re(2048/self.decim)
		self.set_samp_rate(30720e3/self.decim)
		self.set_fft_size(2048/self.decim)

	def get_frames(self):
		return self.frames

	def set_frames(self, frames):
		self.frames = frames

	def get_gain(self):
		return self.gain

	def set_gain(self, gain):
		self.gain = gain

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate
		self.wxgui_scopesink2_0.set_sample_rate(self.samp_rate)
		self.wxgui_scopesink2_0_2.set_sample_rate(self.samp_rate)

	def get_fft_size(self):
		return self.fft_size

	def set_fft_size(self, fft_size):
		self.fft_size = fft_size

	def get_N_re(self):
		return self.N_re

	def set_N_re(self, N_re):
		self.N_re = N_re

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	parser.add_option("", "--N-id-1", dest="N_id_1", type="intx", default=134,
		help="Set N_id_1 [default=%default]")
	parser.add_option("", "--N-id-2", dest="N_id_2", type="intx", default=0,
		help="Set N_id_2 [default=%default]")
	parser.add_option("", "--decim", dest="decim", type="intx", default=16,
		help="Set decim [default=%default]")
	parser.add_option("", "--frames", dest="frames", type="intx", default=1,
		help="Set frames [default=%default]")
	parser.add_option("", "--gain", dest="gain", type="eng_float", default=eng_notation.num_to_str(.08),
		help="Set gain [default=%default]")
	(options, args) = parser.parse_args()
	tb = lte_dl_ss_source(N_id_1=options.N_id_1, N_id_2=options.N_id_2, decim=options.decim, frames=options.frames, gain=options.gain)
	tb.Run(True)

