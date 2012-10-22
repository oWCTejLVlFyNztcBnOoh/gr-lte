#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Pss Corr Block Gui
# Generated: Mon Oct 22 20:34:58 2012
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

	def __init__(self, decim=16, N_id_2=0, avg_frames=1, freq_corr=0, N_id_1=134, file_name="/home/user/git/gr-lte/gr-lte/test/traces/lte_02_796m_30720k_frame.cfile"):
		gr.top_block.__init__(self, "Pss Corr Block Gui")

		##################################################
		# Parameters
		##################################################
		self.decim = decim
		self.N_id_2 = N_id_2
		self.avg_frames = avg_frames
		self.freq_corr = freq_corr
		self.N_id_1 = N_id_1
		self.file_name = file_name

		##################################################
		# Variables
		##################################################
		self.vec_half_frame = vec_half_frame = 30720*5/decim
		self.transition_width = transition_width = 100e3
		self.samp_rate = samp_rate = 30720e3/decim
		self.cutoff_freq = cutoff_freq = 550e3

		##################################################
		# Blocks
		##################################################
		self.gr_vector_sink_sss10 = gr.vector_sink_i(1)
		self.gr_vector_sink_sss0 = gr.vector_sink_i(1)
		self.gr_vector_sink_pss = gr.vector_sink_i(1)
		self.gr_throttle_0 = gr.throttle(gr.sizeof_gr_complex*1, samp_rate*decim)
		self.gr_null_sink_0 = gr.null_sink(gr.sizeof_gr_complex*1)
		self.gr_file_source_0 = gr.file_source(gr.sizeof_gr_complex*1, file_name, False)
		self.fir_filter_xxx_0 = filter.fir_filter_ccc(decim, (firdes.low_pass(1, decim*samp_rate, cutoff_freq, transition_width)))
		self.any_0_0_0 = sss_corr(N_id_1,N_id_2,False, decim, avg_frames,freq_corr)
		self.any_0_0 = sss_corr(N_id_1,N_id_2,True, decim, avg_frames,freq_corr)
		self.any_0 = pss_corr(N_id_2, decim, avg_frames*2,freq_corr)

		##################################################
		# Connections
		##################################################
		self.connect((self.any_0_0_0, 0), (self.gr_vector_sink_sss10, 0))
		self.connect((self.fir_filter_xxx_0, 0), (self.any_0_0_0, 0))
		self.connect((self.any_0_0, 0), (self.gr_vector_sink_sss0, 0))
		self.connect((self.fir_filter_xxx_0, 0), (self.any_0_0, 0))
		self.connect((self.any_0, 0), (self.gr_vector_sink_pss, 0))
		self.connect((self.fir_filter_xxx_0, 0), (self.any_0, 0))
		self.connect((self.gr_throttle_0, 0), (self.fir_filter_xxx_0, 0))
		self.connect((self.gr_file_source_0, 0), (self.gr_throttle_0, 0))
		self.connect((self.gr_file_source_0, 0), (self.gr_null_sink_0, 0))

	def get_decim(self):
		return self.decim

	def set_decim(self, decim):
		self.decim = decim
		self.set_vec_half_frame(30720*5/self.decim)
		self.set_samp_rate(30720e3/self.decim)
		self.fir_filter_xxx_0.set_taps((firdes.low_pass(1, self.decim*self.samp_rate, self.cutoff_freq, self.transition_width)))

	def get_N_id_2(self):
		return self.N_id_2

	def set_N_id_2(self, N_id_2):
		self.N_id_2 = N_id_2

	def get_avg_frames(self):
		return self.avg_frames

	def set_avg_frames(self, avg_frames):
		self.avg_frames = avg_frames

	def get_freq_corr(self):
		return self.freq_corr

	def set_freq_corr(self, freq_corr):
		self.freq_corr = freq_corr

	def get_N_id_1(self):
		return self.N_id_1

	def set_N_id_1(self, N_id_1):
		self.N_id_1 = N_id_1

	def get_file_name(self):
		return self.file_name

	def set_file_name(self, file_name):
		self.file_name = file_name

	def get_vec_half_frame(self):
		return self.vec_half_frame

	def set_vec_half_frame(self, vec_half_frame):
		self.vec_half_frame = vec_half_frame

	def get_transition_width(self):
		return self.transition_width

	def set_transition_width(self, transition_width):
		self.transition_width = transition_width
		self.fir_filter_xxx_0.set_taps((firdes.low_pass(1, self.decim*self.samp_rate, self.cutoff_freq, self.transition_width)))

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate
		self.fir_filter_xxx_0.set_taps((firdes.low_pass(1, self.decim*self.samp_rate, self.cutoff_freq, self.transition_width)))

	def get_cutoff_freq(self):
		return self.cutoff_freq

	def set_cutoff_freq(self, cutoff_freq):
		self.cutoff_freq = cutoff_freq
		self.fir_filter_xxx_0.set_taps((firdes.low_pass(1, self.decim*self.samp_rate, self.cutoff_freq, self.transition_width)))

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	parser.add_option("", "--decim", dest="decim", type="intx", default=16,
		help="Set decim [default=%default]")
	parser.add_option("", "--N-id-2", dest="N_id_2", type="intx", default=0,
		help="Set N_id_2 [default=%default]")
	parser.add_option("", "--avg-frames", dest="avg_frames", type="intx", default=1,
		help="Set avg_frames [default=%default]")
	parser.add_option("", "--freq-corr", dest="freq_corr", type="eng_float", default=eng_notation.num_to_str(0),
		help="Set freq_corr [default=%default]")
	parser.add_option("", "--N-id-1", dest="N_id_1", type="intx", default=134,
		help="Set N_id_1 [default=%default]")
	parser.add_option("", "--file-name", dest="file_name", type="string", default="/home/user/git/gr-lte/gr-lte/test/traces/lte_02_796m_30720k_frame.cfile",
		help="Set file_name [default=%default]")
	(options, args) = parser.parse_args()
	tb = pss_corr_block_gui(decim=options.decim, N_id_2=options.N_id_2, avg_frames=options.avg_frames, freq_corr=options.freq_corr, N_id_1=options.N_id_1, file_name=options.file_name)
	tb.run()

