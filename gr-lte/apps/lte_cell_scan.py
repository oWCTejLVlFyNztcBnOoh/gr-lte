#!/usr/bin/env python

from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.gr import firdes
from optparse import OptionParser
from pss_corr import *  
from sss_corr import *
from lte_dl_ss_source import *   
import logging

class lte_cell_scan:
	"""
	Cell scanner for LTE signals
	
	concatenate((zeros((160+144*4+2048*5)/decim), gen_sss_td(N_id_1, N_id_2, True, N_re).get_data(), 
	  gen_pss_td(N_id_2, N_re, 144).get_data(), zeros((30720/2 + 30720*4)/decim), 
	  zeros((160+144*4+2048*5)/decim), gen_sss_td(N_id_1, N_id_2, False, N_re).get_data(), 
	  gen_pss_td(N_id_2, N_re, 144).get_data(), zeros((30720/2 + 30720*4)/decim)),1)
	  
	/home/user/git/gr-lte/gr-lte/test/traces/lte_02_796m_30720k_frame.cfile
	"""
	
	def __init__(self, file_name, decim=16, avg_frames=1, freq_corr=0, dump=None):
		##################################################
		# Parameters
		##################################################
		self.decim = decim
		self.avg_frames = avg_frames
		self.freq_corr = freq_corr
		self.file_name = file_name
		self.dump = dump

		##################################################
		# Variables
		##################################################
		self.vec_half_frame = vec_half_frame = 30720*5/decim
		self.samp_rate = samp_rate = 30720e3/decim
		self.cutoff_freq = cutoff_freq = 550e3
		self.transition_width = transition_width = 100e3

		##################################################
		# Generate input vector source
		##################################################
		if None != file_name:
			logging.debug("Reading from file: " + file_name)
			file_source = gr.file_source(gr.sizeof_gr_complex*1, file_name, False)
		else:
			N_id_1=134
			N_id_2=0
			logging.debug("Generating DL sync signals for cell ID {}".format(3*N_id_1+N_id_2))
			file_source = lte_dl_ss_source(decim=1,N_id_1=N_id_1,N_id_2=N_id_2,frames=self.avg_frames)
		decim_lowpass = filter.fir_filter_ccc(decim, (firdes.low_pass(1, decim*samp_rate, cutoff_freq, transition_width)))
		sink = gr.vector_sink_c();
		top = gr.top_block("input reader graph")
		top.connect(file_source, decim_lowpass, sink)
		if dump != None:
			file_sink = gr.file_sink(gr.sizeof_gr_complex, self.dump + "_input.cfile")
			top.connect(decim_lowpass, file_sink)
		top.run()
		logging.debug("No. samples to process: {}".format(len(sink.data()))) 
		
		if -1 == self.avg_frames:
			self.avg_frames = int(floor(len(sink.data()) / vec_half_frame / 2))
		logging.debug("Config: decim {}, avg_frames {}, freq_corr {}".format(self.decim, self.avg_frames, self.freq_corr)) 
		self.source = gr.vector_source_c(sink.data())
		

	def scan(self):
		res = self.correlate_pss()
		logging.debug("Found PSS correlations " + str(res))
		for (N_id_2, frame_time, peak) in res:
			self.correlate_sss(N_id_2, frame_time, peak)


	def correlate_pss(self):
		top = gr.top_block("pss corr graph")
		#and_block = gr.add_cc()
		#top.connect(self.source, and_block)
		corr = range(0,3)
		sink_pss = range(0,3)
		for N_id_2 in range(0,3):
			corr[N_id_2] = pss_corr(N_id_2, self.decim, self.avg_frames*2, self.freq_corr, self.dump)
			sink_pss[N_id_2] = gr.vector_sink_f()
			top.connect(self.source, corr[N_id_2], sink_pss[N_id_2])
		top.run()
		ret = []
		for N_id_2 in range(0,3):
			pss_res = sink_pss[N_id_2].data()
			for i in range(0, len(pss_res)/2): 
				ret += [(N_id_2, pss_res[2*i], pss_res[2*i+1])]
		return ret
			

	def correlate_sss(self, N_id_2, frame_time, pss_peak):
		# TODO ofdm_receiver.py
		logging.debug("Searching for SSS: N_id_2 {}, frame time {:5.0f}, corr peak {:7.0f}".format(N_id_2, frame_time, pss_peak))
		#self.gr_vector_sink_sss10 = gr.vector_sink_i(1)
		#self.gr_vector_sink_sss0 = gr.vector_sink_i(1)
		#self.any_0_0_0 = sss_corr(N_id_1,N_id_2,False, decim, avg_frames,freq_corr)
		#self.any_0_0 = sss_corr(N_id_1,N_id_2,True, decim, avg_frames,freq_corr)
		#if logging:
		#	self.connect(self.chan_filt, gr.file_sink(gr.sizeof_gr_complex, "ofdm_receiver-chan_filt_c.dat"))

		pass


if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	parser.add_option("", "--decim", dest="decim", type="intx", default=16,
		help="Set decim [default=%default]")
	parser.add_option("", "--avg-frames", dest="avg_frames", type="intx", default=-1,
		help="Set avg_frames [default=%default]")
	parser.add_option("", "--freq-corr", dest="freq_corr", type="eng_float", default=eng_notation.num_to_str(0),
		help="Set freq_corr [default=%default]")
	parser.add_option("", "--file-name", dest="file_name", type="string", default=None,
		help="Set file_name [default=%default]")
	parser.add_option("", "--dump", dest="dump", type="string", default=None,
		help="dump intermediate results [default=%default]")
	(options, args) = parser.parse_args()
	logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)
	tb = lte_cell_scan(decim=options.decim, avg_frames=options.avg_frames, freq_corr=options.freq_corr, file_name=options.file_name, dump=options.dump)
	tb.scan()


					
					