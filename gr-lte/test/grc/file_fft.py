#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: File Fft
# Generated: Sat Oct 27 16:54:13 2012
##################################################

from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import window
from gnuradio.eng_option import eng_option
from gnuradio.gr import firdes
from gnuradio.wxgui import scopesink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
import wx

class file_fft(grc_wxgui.top_block_gui):

	def __init__(self, decim=16):
		grc_wxgui.top_block_gui.__init__(self, title="File Fft")
		_icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
		self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

		##################################################
		# Parameters
		##################################################
		self.decim = decim

		##################################################
		# Variables
		##################################################
		self.samp_rate = samp_rate = 32000
		self.fft_size = fft_size = 2048/decim
		self.N_re = N_re = 62

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
		self.gr_vector_to_stream_0 = gr.vector_to_stream(gr.sizeof_gr_complex*1, fft_size)
		self.gr_stream_to_vector_0 = gr.stream_to_vector(gr.sizeof_gr_complex*1, fft_size)
		self.gr_keep_m_in_n_0 = gr.keep_m_in_n(gr.sizeof_gr_complex, N_re, fft_size, (fft_size-N_re)/2)
		self.gr_file_source_0 = gr.file_source(gr.sizeof_gr_complex*1, "/home/user/git/gr-lte/gr-lte/test/foo_pss0_sss_in.cfile", True)
		self.gr_fft_vxx_0 = gr.fft_vcc(fft_size, True, (window.blackmanharris(1024)), True, 1)

		##################################################
		# Connections
		##################################################
		self.connect((self.gr_stream_to_vector_0, 0), (self.gr_fft_vxx_0, 0))
		self.connect((self.gr_file_source_0, 0), (self.gr_stream_to_vector_0, 0))
		self.connect((self.gr_vector_to_stream_0, 0), (self.gr_keep_m_in_n_0, 0))
		self.connect((self.gr_fft_vxx_0, 0), (self.gr_vector_to_stream_0, 0))
		self.connect((self.gr_keep_m_in_n_0, 0), (self.wxgui_scopesink2_0, 0))

	def get_decim(self):
		return self.decim

	def set_decim(self, decim):
		self.decim = decim
		self.set_fft_size(2048/self.decim)

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate
		self.wxgui_scopesink2_0.set_sample_rate(self.samp_rate)

	def get_fft_size(self):
		return self.fft_size

	def set_fft_size(self, fft_size):
		self.fft_size = fft_size
		self.gr_keep_m_in_n_0.set_offset((self.fft_size-self.N_re)/2)
		self.gr_keep_m_in_n_0.set_n(self.fft_size)

	def get_N_re(self):
		return self.N_re

	def set_N_re(self, N_re):
		self.N_re = N_re
		self.gr_keep_m_in_n_0.set_offset((self.fft_size-self.N_re)/2)
		self.gr_keep_m_in_n_0.set_m(self.N_re)

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	parser.add_option("", "--decim", dest="decim", type="intx", default=16,
		help="Set decim [default=%default]")
	(options, args) = parser.parse_args()
	tb = file_fft(decim=options.decim)
	tb.Run(True)

