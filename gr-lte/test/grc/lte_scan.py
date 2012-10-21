#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Lte Scan
# Generated: Thu Oct 18 20:22:38 2012
##################################################

from gnuradio import LTE_fdd_dl_fs
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.gr import firdes
from optparse import OptionParser

class lte_scan(gr.top_block):

	def __init__(self, center_freq=-12421):
		gr.top_block.__init__(self, "Lte Scan")

		##################################################
		# Parameters
		##################################################
		self.center_freq = center_freq

		##################################################
		# Variables
		##################################################
		self.samp_rate = samp_rate = 2160e3
		self.lte_rate = lte_rate = 30.72e6

		##################################################
		# Blocks
		##################################################
		self.gr_interleave_0 = gr.interleave(gr.sizeof_float*1)
		self.gr_float_to_char_0 = gr.float_to_char(1, 1)
		self.gr_file_source_0 = gr.file_source(gr.sizeof_gr_complex*1, "/home/user/git/mbus/gsm-test/capture/lte_02_796m_30720k_frame.cfile", False)
		self.gr_complex_to_float_0 = gr.complex_to_float(1)
		self.freq_xlating_fir_filter_xxx_0 = filter.freq_xlating_fir_filter_ccc(1, (1, ), center_freq, lte_rate)
		self.any_sink_0 = LTE_fdd_dl_fs.samp_buf()

		##################################################
		# Connections
		##################################################
		self.connect((self.gr_complex_to_float_0, 0), (self.gr_interleave_0, 0))
		self.connect((self.gr_complex_to_float_0, 1), (self.gr_interleave_0, 1))
		self.connect((self.gr_interleave_0, 0), (self.gr_float_to_char_0, 0))
		self.connect((self.gr_float_to_char_0, 0), (self.any_sink_0, 0))
		self.connect((self.freq_xlating_fir_filter_xxx_0, 0), (self.gr_complex_to_float_0, 0))
		self.connect((self.gr_file_source_0, 0), (self.freq_xlating_fir_filter_xxx_0, 0))

	def get_center_freq(self):
		return self.center_freq

	def set_center_freq(self, center_freq):
		self.center_freq = center_freq
		self.freq_xlating_fir_filter_xxx_0.set_center_freq(self.center_freq)

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate

	def get_lte_rate(self):
		return self.lte_rate

	def set_lte_rate(self, lte_rate):
		self.lte_rate = lte_rate

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	parser.add_option("-c", "--center-freq", dest="center_freq", type="eng_float", default=eng_notation.num_to_str(-12421),
		help="Set center_freq [default=%default]")
	(options, args) = parser.parse_args()
	tb = lte_scan(center_freq=options.center_freq)
	tb.run()

