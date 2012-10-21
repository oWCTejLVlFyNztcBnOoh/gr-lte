#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Lte Cp
# Generated: Sun Oct 14 18:50:10 2012
##################################################

from baz import facsink
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.gr import firdes
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
import wx

class lte_cp(grc_wxgui.top_block_gui):

	def __init__(self):
		grc_wxgui.top_block_gui.__init__(self, title="Lte Cp")
		_icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
		self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

		##################################################
		# Variables
		##################################################
		self.samp_rate = samp_rate = 2160e3
		self.lte_rate = lte_rate = 30.72e6

		##################################################
		# Blocks
		##################################################
		self.gr_throttle_0 = gr.throttle(gr.sizeof_gr_complex*1, samp_rate)
		self.gr_fractional_interpolator_xx_0 = gr.fractional_interpolator_cc(0, samp_rate/lte_rate)
		self.gr_file_source_0 = gr.file_source(gr.sizeof_gr_complex*1, "/home/user/git/mbus/gsm-test/capture/lte_02_796m_2160k.cfile", True)
		self.facsink_0 = facsink.fac_sink_c(
			self.GetWin(),
			title="Fast AutoCorrelation",
			sample_rate=lte_rate,
			baseband_freq=0,
		        y_per_div=10,
			ref_level=50,
			fac_size=512*16,
		        fac_rate=facsink.default_fac_rate,
		        average=False,
			avg_alpha=0,
			peak_hold=False,
		)
		self.Add(self.facsink_0.win)

		##################################################
		# Connections
		##################################################
		self.connect((self.gr_file_source_0, 0), (self.gr_throttle_0, 0))
		self.connect((self.gr_throttle_0, 0), (self.gr_fractional_interpolator_xx_0, 0))
		self.connect((self.gr_fractional_interpolator_xx_0, 0), (self.facsink_0, 0))

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate
		self.gr_fractional_interpolator_xx_0.set_interp_ratio(self.samp_rate/self.lte_rate)

	def get_lte_rate(self):
		return self.lte_rate

	def set_lte_rate(self, lte_rate):
		self.lte_rate = lte_rate
		self.gr_fractional_interpolator_xx_0.set_interp_ratio(self.samp_rate/self.lte_rate)
		self.facsink_0.set_sample_rate(self.lte_rate)

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	(options, args) = parser.parse_args()
	tb = lte_cp()
	tb.Run(True)

