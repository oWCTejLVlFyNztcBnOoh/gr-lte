#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Gen Pss Gui
# Generated: Sun Oct 21 21:34:45 2012
##################################################

from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import window
from gnuradio.eng_option import eng_option
from gnuradio.gr import firdes
from gnuradio.wxgui import fftsink2
from gnuradio.wxgui import scopesink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
from pss_source import *
from sss_source import * 
import wx

class gen_pss_gui(grc_wxgui.top_block_gui):

	def __init__(self, N_cp_ts=144, freq_corr=0, N_re=128):
		grc_wxgui.top_block_gui.__init__(self, title="Gen Pss Gui")
		_icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
		self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

		##################################################
		# Parameters
		##################################################
		self.N_cp_ts = N_cp_ts
		self.freq_corr = freq_corr
		self.N_re = N_re

		##################################################
		# Variables
		##################################################
		self.samp_rate = samp_rate = 30.72e6/2048*N_re

		##################################################
		# Blocks
		##################################################
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
		self.wxgui_fftsink2_0 = fftsink2.fft_sink_c(
			self.GetWin(),
			baseband_freq=0,
			y_per_div=10,
			y_divs=10,
			ref_level=0,
			ref_scale=2.0,
			sample_rate=samp_rate,
			fft_size=1024,
			fft_rate=15,
			average=False,
			avg_alpha=None,
			title="FFT Plot",
			peak_hold=False,
		)
		self.Add(self.wxgui_fftsink2_0.win)
		self.gr_vector_source_x_0 = gr.vector_source_c((gen_sss_td(0, 134, False).get_data()), True, 1)
		self.gr_throttle_0 = gr.throttle(gr.sizeof_gr_complex*1, samp_rate)

		##################################################
		# Connections
		##################################################
		self.connect((self.gr_throttle_0, 0), (self.wxgui_scopesink2_0, 0))
		self.connect((self.gr_throttle_0, 0), (self.wxgui_fftsink2_0, 0))
		self.connect((self.gr_vector_source_x_0, 0), (self.gr_throttle_0, 0))

	def get_N_cp_ts(self):
		return self.N_cp_ts

	def set_N_cp_ts(self, N_cp_ts):
		self.N_cp_ts = N_cp_ts

	def get_freq_corr(self):
		return self.freq_corr

	def set_freq_corr(self, freq_corr):
		self.freq_corr = freq_corr

	def get_N_re(self):
		return self.N_re

	def set_N_re(self, N_re):
		self.N_re = N_re
		self.set_samp_rate(30.72e6/2048*self.N_re)

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate
		self.wxgui_fftsink2_0.set_sample_rate(self.samp_rate)
		self.wxgui_scopesink2_0.set_sample_rate(self.samp_rate)

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	parser.add_option("", "--N-cp-ts", dest="N_cp_ts", type="intx", default=144,
		help="Set N_cp_ts [default=%default]")
	parser.add_option("", "--freq-corr", dest="freq_corr", type="eng_float", default=eng_notation.num_to_str(0),
		help="Set freq_corr [default=%default]")
	parser.add_option("", "--N-re", dest="N_re", type="intx", default=128,
		help="Set N_re [default=%default]")
	(options, args) = parser.parse_args()
	tb = gen_pss_gui(N_cp_ts=options.N_cp_ts, freq_corr=options.freq_corr, N_re=options.N_re)
	tb.Run(True)

