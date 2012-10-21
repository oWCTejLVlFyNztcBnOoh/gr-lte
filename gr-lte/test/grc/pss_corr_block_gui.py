#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Pss Corr Block Gui
# Generated: Sun Oct 21 18:38:01 2012
##################################################

from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.gr import firdes
from optparse import OptionParser
from pss_corr import *  

class pss_corr_block_gui(gr.top_block):

	def __init__(self, N_id_2=0, decim=16, avg_hf=1):
		gr.top_block.__init__(self, "Pss Corr Block Gui")

		##################################################
		# Parameters
		##################################################
		self.N_id_2 = N_id_2
		self.decim = decim
		self.avg_hf = avg_hf

		##################################################
		# Variables
		##################################################
		self.vec_half_frame = vec_half_frame = 30720*5/decim
		self.samp_rate = samp_rate = 30720e3/decim

		##################################################
		# Blocks
		##################################################
		self.gr_vector_sink_x_0 = gr.vector_sink_i(1)
		self.gr_throttle_0 = gr.throttle(gr.sizeof_gr_complex*1, samp_rate*decim)
		self.gr_file_source_0 = gr.file_source(gr.sizeof_gr_complex*1, "/home/user/git/gr-lte/gr-lte/test/traces/lte_02_796m_30720k_frame.cfile", False)
		self.fir_filter_xxx_0 = filter.fir_filter_ccc(decim, (firdes.low_pass(1, decim*samp_rate, 550e3, 100e3)))
		self.any_0 = pss_corr(N_id_2, decim, avg_hf)

		##################################################
		# Connections
		##################################################
		self.connect((self.gr_file_source_0, 0), (self.gr_throttle_0, 0))
		self.connect((self.gr_throttle_0, 0), (self.fir_filter_xxx_0, 0))
		self.connect((self.fir_filter_xxx_0, 0), (self.any_0, 0))
		self.connect((self.any_0, 0), (self.gr_vector_sink_x_0, 0))

	def get_N_id_2(self):
		return self.N_id_2

	def set_N_id_2(self, N_id_2):
		self.N_id_2 = N_id_2

	def get_decim(self):
		return self.decim

	def set_decim(self, decim):
		self.decim = decim
		self.fir_filter_xxx_0.set_taps((firdes.low_pass(1, self.decim*self.samp_rate, 550e3, 100e3)))
		self.set_vec_half_frame(30720*5/self.decim)
		self.set_samp_rate(30720e3/self.decim)

	def get_avg_hf(self):
		return self.avg_hf

	def set_avg_hf(self, avg_hf):
		self.avg_hf = avg_hf

	def get_vec_half_frame(self):
		return self.vec_half_frame

	def set_vec_half_frame(self, vec_half_frame):
		self.vec_half_frame = vec_half_frame

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate
		self.fir_filter_xxx_0.set_taps((firdes.low_pass(1, self.decim*self.samp_rate, 550e3, 100e3)))

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	parser.add_option("", "--N-id-2", dest="N_id_2", type="intx", default=0,
		help="Set N_id_2 [default=%default]")
	parser.add_option("", "--decim", dest="decim", type="intx", default=16,
		help="Set decim [default=%default]")
	parser.add_option("", "--avg-hf", dest="avg_hf", type="intx", default=1,
		help="Set avg_hf [default=%default]")
	(options, args) = parser.parse_args()
	tb = pss_corr_block_gui(N_id_2=options.N_id_2, decim=options.decim, avg_hf=options.avg_hf)
	tb.run()
	for f in tb.gr_vector_sink_x_0.data():
		print f

