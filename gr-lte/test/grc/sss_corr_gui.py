#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Sss Corr Gui
# Generated: Mon Oct 22 21:52:21 2012
##################################################

from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio import window
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.gr import firdes
from gnuradio.wxgui import scopesink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
from sss_source import *
import wx

class sss_corr_gui(grc_wxgui.top_block_gui):

	def __init__(self, decim=16, N_id_1=134, N_id_2=0):
		grc_wxgui.top_block_gui.__init__(self, title="Sss Corr Gui")
		_icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
		self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

		##################################################
		# Parameters
		##################################################
		self.decim = decim
		self.N_id_1 = N_id_1
		self.N_id_2 = N_id_2

		##################################################
		# Variables
		##################################################
		self.sss_start_ts = sss_start_ts = 10608 - 2048 - 144
		self.samp_rate = samp_rate = 30720e3/decim

		##################################################
		# Blocks
		##################################################
		self.wxgui_scopesink2_0 = scopesink2.scope_sink_f(
			self.GetWin(),
			title="Scope Plot",
			sample_rate=200,
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
		self.gr_vector_to_stream_0 = gr.vector_to_stream(gr.sizeof_gr_complex*1, 2048/decim)
		self.gr_throttle_0 = gr.throttle(gr.sizeof_gr_complex*1, samp_rate*decim)
		self.gr_stream_to_vector_0 = gr.stream_to_vector(gr.sizeof_gr_complex*1, 2048/decim)
		self.gr_skiphead_0 = gr.skiphead(gr.sizeof_gr_complex*1, sss_start_ts)
		self.gr_multiply_const_vxx_0 = gr.multiply_const_vcc((gen_sss_fd(N_id_1, N_id_2, 2048/decim).get_sss_conj_rev()))
		self.gr_integrate_xx_0 = gr.integrate_cc(2048/decim)
		self.gr_file_source_0 = gr.file_source(gr.sizeof_gr_complex*1, "/home/user/git/gr-lte/gr-lte/test/traces/lte_02_796m_30720k_frame.cfile", True)
		self.gr_fft_vxx_0 = gr.fft_vcc(2048/decim, True, (window.blackmanharris(1024)), True, 1)
		self.gr_complex_to_mag_0 = gr.complex_to_mag(1)
		self.fir_filter_xxx_0 = filter.fir_filter_ccc(decim, (firdes.low_pass(1, decim*samp_rate, 550e3, 100e3)))

		##################################################
		# Connections
		##################################################
		self.connect((self.gr_file_source_0, 0), (self.gr_throttle_0, 0))
		self.connect((self.gr_stream_to_vector_0, 0), (self.gr_fft_vxx_0, 0))
		self.connect((self.gr_fft_vxx_0, 0), (self.gr_multiply_const_vxx_0, 0))
		self.connect((self.gr_multiply_const_vxx_0, 0), (self.gr_vector_to_stream_0, 0))
		self.connect((self.gr_vector_to_stream_0, 0), (self.gr_integrate_xx_0, 0))
		self.connect((self.gr_integrate_xx_0, 0), (self.gr_complex_to_mag_0, 0))
		self.connect((self.gr_complex_to_mag_0, 0), (self.wxgui_scopesink2_0, 0))
		self.connect((self.gr_throttle_0, 0), (self.gr_skiphead_0, 0))
		self.connect((self.gr_skiphead_0, 0), (self.fir_filter_xxx_0, 0))
		self.connect((self.fir_filter_xxx_0, 0), (self.gr_stream_to_vector_0, 0))

	def get_decim(self):
		return self.decim

	def set_decim(self, decim):
		self.decim = decim
		self.set_samp_rate(30720e3/self.decim)
		self.fir_filter_xxx_0.set_taps((firdes.low_pass(1, self.decim*self.samp_rate, 550e3, 100e3)))
		self.gr_multiply_const_vxx_0.set_k((gen_sss_fd(self.N_id_1, self.N_id_2, 2048/self.decim).get_sss_conj_rev()))

	def get_N_id_1(self):
		return self.N_id_1

	def set_N_id_1(self, N_id_1):
		self.N_id_1 = N_id_1
		self.gr_multiply_const_vxx_0.set_k((gen_sss_fd(self.N_id_1, self.N_id_2, 2048/self.decim).get_sss_conj_rev()))

	def get_N_id_2(self):
		return self.N_id_2

	def set_N_id_2(self, N_id_2):
		self.N_id_2 = N_id_2
		self.gr_multiply_const_vxx_0.set_k((gen_sss_fd(self.N_id_1, self.N_id_2, 2048/self.decim).get_sss_conj_rev()))

	def get_sss_start_ts(self):
		return self.sss_start_ts

	def set_sss_start_ts(self, sss_start_ts):
		self.sss_start_ts = sss_start_ts

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate
		self.fir_filter_xxx_0.set_taps((firdes.low_pass(1, self.decim*self.samp_rate, 550e3, 100e3)))

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	parser.add_option("", "--N-id-1", dest="N_id_1", type="intx", default=134,
		help="Set N_id_1 [default=%default]")
	parser.add_option("", "--N-id-2", dest="N_id_2", type="intx", default=0,
		help="Set N_id_2 [default=%default]")
	(options, args) = parser.parse_args()
	tb = sss_corr_gui(N_id_1=options.N_id_1, N_id_2=options.N_id_2)
	tb.Run(True)

