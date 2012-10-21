#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Pss Corr Gui
# Generated: Sun Oct 21 17:43:05 2012
##################################################

from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.gr import firdes
from gnuradio.wxgui import numbersink2
from gnuradio.wxgui import scopesink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
from pss_source import *
import baz
import wx

class pss_corr_gui(grc_wxgui.top_block_gui):

	def __init__(self, decim=16, N_id_2=2):
		grc_wxgui.top_block_gui.__init__(self, title="Pss Corr Gui")
		_icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
		self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

		##################################################
		# Parameters
		##################################################
		self.decim = decim
		self.N_id_2 = N_id_2

		##################################################
		# Variables
		##################################################
		self.vec_half_frame = vec_half_frame = 30720*5/decim
		self.samp_rate = samp_rate = 30720e3/decim

		##################################################
		# Blocks
		##################################################
		self.wxgui_scopesink2_0_0_0 = scopesink2.scope_sink_f(
			self.GetWin(),
			title="Scope Plot",
			sample_rate=samp_rate/vec_half_frame,
			v_scale=0,
			v_offset=0,
			t_scale=0,
			ac_couple=False,
			xy_mode=False,
			num_inputs=1,
			trig_mode=gr.gr_TRIG_MODE_AUTO,
			y_axis_label="Counts",
		)
		self.Add(self.wxgui_scopesink2_0_0_0.win)
		self.wxgui_numbersink2_0 = numbersink2.number_sink_f(
			self.GetWin(),
			unit="Units",
			minval=0,
			maxval=vec_half_frame*decim,
			factor=1.0,
			decimal_places=0,
			ref_level=0,
			sample_rate=samp_rate/vec_half_frame,
			number_rate=15,
			average=False,
			avg_alpha=None,
			label="Number Plot",
			peak_hold=False,
			show_gauge=True,
		)
		self.Add(self.wxgui_numbersink2_0.win)
		self.gr_throttle_0 = gr.throttle(gr.sizeof_gr_complex*1, samp_rate*decim)
		self.gr_stream_to_vector_0 = gr.stream_to_vector(gr.sizeof_float*1, vec_half_frame)
		self.gr_short_to_float_0 = gr.short_to_float(1, 1)
		self.gr_null_sink_0 = gr.null_sink(gr.sizeof_short*1)
		self.gr_multiply_const_vxx_0 = gr.multiply_const_vff((decim, ))
		self.gr_file_source_0 = gr.file_source(gr.sizeof_gr_complex*1, "/home/user/git/gr-lte/gr-lte/test/traces/lte_02_796m_30720k_frame.cfile", True)
		self.gr_deinterleave_0 = gr.deinterleave(gr.sizeof_float*vec_half_frame)
		self.gr_complex_to_mag_0 = gr.complex_to_mag(1)
		self.gr_argmax_xx_0 = gr.argmax_fs(vec_half_frame)
		self.gr_add_xx_0 = gr.add_vff(vec_half_frame)
		self.fir_filter_xxx_0_0 = filter.fir_filter_ccc(1, (gen_pss_td(N_id_2).get_data_conj_rev()))
		self.fir_filter_xxx_0 = filter.fir_filter_ccc(decim, (firdes.low_pass(1, decim*samp_rate, 550e3, 100e3)))

		##################################################
		# Connections
		##################################################
		self.connect((self.fir_filter_xxx_0, 0), (self.fir_filter_xxx_0_0, 0))
		self.connect((self.fir_filter_xxx_0_0, 0), (self.gr_complex_to_mag_0, 0))
		self.connect((self.gr_file_source_0, 0), (self.gr_throttle_0, 0))
		self.connect((self.gr_throttle_0, 0), (self.fir_filter_xxx_0, 0))
		self.connect((self.gr_complex_to_mag_0, 0), (self.gr_stream_to_vector_0, 0))
		self.connect((self.gr_stream_to_vector_0, 0), (self.gr_deinterleave_0, 0))
		self.connect((self.gr_deinterleave_0, 0), (self.gr_add_xx_0, 0))
		self.connect((self.gr_deinterleave_0, 1), (self.gr_add_xx_0, 1))
		self.connect((self.gr_deinterleave_0, 2), (self.gr_add_xx_0, 2))
		self.connect((self.gr_deinterleave_0, 3), (self.gr_add_xx_0, 3))
		self.connect((self.gr_argmax_xx_0, 1), (self.gr_null_sink_0, 0))
		self.connect((self.gr_argmax_xx_0, 0), (self.gr_short_to_float_0, 0))
		self.connect((self.gr_short_to_float_0, 0), (self.gr_multiply_const_vxx_0, 0))
		self.connect((self.gr_multiply_const_vxx_0, 0), (self.wxgui_numbersink2_0, 0))
		self.connect((self.gr_multiply_const_vxx_0, 0), (self.wxgui_scopesink2_0_0_0, 0))
		self.connect((self.gr_add_xx_0, 0), (self.gr_argmax_xx_0, 0))

	def get_decim(self):
		return self.decim

	def set_decim(self, decim):
		self.decim = decim
		self.fir_filter_xxx_0.set_taps((firdes.low_pass(1, self.decim*self.samp_rate, 550e3, 100e3)))
		self.set_vec_half_frame(30720*5/self.decim)
		self.set_samp_rate(30720e3/self.decim)
		self.gr_multiply_const_vxx_0.set_k((self.decim, ))

	def get_N_id_2(self):
		return self.N_id_2

	def set_N_id_2(self, N_id_2):
		self.N_id_2 = N_id_2
		self.fir_filter_xxx_0_0.set_taps((gen_pss_td(self.N_id_2).get_data_conj_rev()))

	def get_vec_half_frame(self):
		return self.vec_half_frame

	def set_vec_half_frame(self, vec_half_frame):
		self.vec_half_frame = vec_half_frame
		self.wxgui_scopesink2_0_0_0.set_sample_rate(self.samp_rate/self.vec_half_frame)

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate
		self.fir_filter_xxx_0.set_taps((firdes.low_pass(1, self.decim*self.samp_rate, 550e3, 100e3)))
		self.wxgui_scopesink2_0_0_0.set_sample_rate(self.samp_rate/self.vec_half_frame)

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	parser.add_option("", "--N-id-2", dest="N_id_2", type="intx", default=2,
		help="Set N_id_2 [default=%default]")
	(options, args) = parser.parse_args()
	tb = pss_corr_gui(N_id_2=options.N_id_2)
	tb.Run(True)

