#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Pss Corr Block Gui
# Generated: Sun Oct 21 21:54:33 2012
##################################################

from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.gr import firdes
from optparse import OptionParser
from pss_corr import *  
from sss_corr import *   

class pss_corr_block_gui(gr.top_block):

	def __init__(self, decim=16, N_id_2=0, N_id_1=134, avg_frames=1):
		gr.top_block.__init__(self, "Pss Corr Block Gui")

		##################################################
		# Parameters
		##################################################
		self.decim = decim
		self.N_id_2 = N_id_2
		self.N_id_1 = N_id_1
		self.avg_frames = avg_frames

		##################################################
		# Variables
		##################################################
		self.vec_half_frame = vec_half_frame = 30720*5/decim
		self.samp_rate = samp_rate = 30720e3/decim

		##################################################
		# Blocks
		##################################################
		self.gr_vector_sink_sss10 = gr.vector_sink_i(1)
		self.gr_vector_sink_sss0 = gr.vector_sink_i(1)
		self.gr_vector_sink_pss = gr.vector_sink_i(1)
		self.gr_throttle_0 = gr.throttle(gr.sizeof_gr_complex*1, samp_rate*decim)
		self.gr_file_source_0 = gr.file_source(gr.sizeof_gr_complex*1, "/home/user/git/gr-lte/gr-lte/test/traces/lte_02_796m_30720k_frame.cfile", False)
		self.fir_filter_xxx_0 = filter.fir_filter_ccc(decim, (firdes.low_pass(1, decim*samp_rate, 550e3, 100e3)))
		self.any_0_0_0 = sss_corr(N_id_1,N_id_2,False, decim, avg_frames)
		self.any_0_0 = sss_corr(N_id_1,N_id_2,True, decim, avg_frames)
		self.any_0 = pss_corr(N_id_2, decim, avg_frames*2)

		##################################################
		# Connections
		##################################################
		self.connect((self.gr_file_source_0, 0), (self.gr_throttle_0, 0))
		self.connect((self.gr_throttle_0, 0), (self.fir_filter_xxx_0, 0))
		self.connect((self.fir_filter_xxx_0, 0), (self.any_0, 0))
		self.connect((self.any_0, 0), (self.gr_vector_sink_pss, 0))
		self.connect((self.fir_filter_xxx_0, 0), (self.any_0_0, 0))
		self.connect((self.any_0_0, 0), (self.gr_vector_sink_sss0, 0))
		self.connect((self.fir_filter_xxx_0, 0), (self.any_0_0_0, 0))
		self.connect((self.any_0_0_0, 0), (self.gr_vector_sink_sss10, 0))

	def get_decim(self):
		return self.decim

	def set_decim(self, decim):
		self.decim = decim
		self.set_vec_half_frame(30720*5/self.decim)
		self.set_samp_rate(30720e3/self.decim)
		self.fir_filter_xxx_0.set_taps((firdes.low_pass(1, self.decim*self.samp_rate, 550e3, 100e3)))

	def get_N_id_2(self):
		return self.N_id_2

	def set_N_id_2(self, N_id_2):
		self.N_id_2 = N_id_2

	def get_N_id_1(self):
		return self.N_id_1

	def set_N_id_1(self, N_id_1):
		self.N_id_1 = N_id_1

	def get_avg_frames(self):
		return self.avg_frames

	def set_avg_frames(self, avg_frames):
		self.avg_frames = avg_frames

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
	parser.add_option("", "--decim", dest="decim", type="intx", default=16,
		help="Set decim [default=%default]")
	parser.add_option("", "--N-id-2", dest="N_id_2", type="intx", default=0,
		help="Set N_id_2 [default=%default]")
	parser.add_option("", "--N-id-1", dest="N_id_1", type="intx", default=134,
		help="Set N_id_1 [default=%default]")
	parser.add_option("", "--avg-frames", dest="avg_frames", type="intx", default=1,
		help="Set avg_frames [default=%default]")
	(options, args) = parser.parse_args()
	tb = pss_corr_block_gui(decim=options.decim, N_id_2=options.N_id_2, N_id_1=options.N_id_1, avg_frames=options.avg_frames)
	tb.run()
	print_result(tb)
	#pss = tb.gr_vector_sink_pss.data()
	#sss0 = tb.gr_vector_sink_sss0.data()
	#sss10 = tb.gr_vector_sink_sss10.data()
	#for i in range(0,len(pss)):
	#  print "{:06d} {:06d} {:06d} {:06d}".format(pss[i], pss[i]+30720*5, sss0[i], sss10[i]) 
    