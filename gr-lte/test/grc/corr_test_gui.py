#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Corr Test Gui
# Generated: Mon Oct 22 21:06:32 2012
##################################################

from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.gr import firdes
from gnuradio.wxgui import scopesink2
from grc_gnuradio import wxgui as grc_wxgui
from numpy import *
from optparse import OptionParser
import wx

class corr_test_gui(grc_wxgui.top_block_gui):

	def __init__(self):
		grc_wxgui.top_block_gui.__init__(self, title="Corr Test Gui")
		_icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
		self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

		##################################################
		# Variables
		##################################################
		self.samp_rate = samp_rate = 32000

		##################################################
		# Blocks
		##################################################
		self.wxgui_scopesink2_0 = scopesink2.scope_sink_f(
			self.GetWin(),
			title="Scope Plot",
			sample_rate=samp_rate,
			v_scale=0,
			v_offset=0,
			t_scale=0,
			ac_couple=False,
			xy_mode=False,
			num_inputs=2,
			trig_mode=gr.gr_TRIG_MODE_AUTO,
			y_axis_label="Counts",
		)
		self.Add(self.wxgui_scopesink2_0.win)
		self.gr_vector_source_x_0 = gr.vector_source_f((concatenate((ones(10),zeros(200)),1)), True, 1)
		self.gr_throttle_0 = gr.throttle(gr.sizeof_float*1, samp_rate)
		self.fir_filter_xxx_0 = filter.fir_filter_fff(1, (ones(10)))

		##################################################
		# Connections
		##################################################
		self.connect((self.gr_vector_source_x_0, 0), (self.gr_throttle_0, 0))
		self.connect((self.gr_throttle_0, 0), (self.fir_filter_xxx_0, 0))
		self.connect((self.fir_filter_xxx_0, 0), (self.wxgui_scopesink2_0, 0))
		self.connect((self.gr_throttle_0, 0), (self.wxgui_scopesink2_0, 1))

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate
		self.wxgui_scopesink2_0.set_sample_rate(self.samp_rate)

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	(options, args) = parser.parse_args()
	tb = corr_test_gui()
	tb.Run(True)

