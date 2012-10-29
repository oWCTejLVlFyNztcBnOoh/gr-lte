#!/usr/bin/env python

from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.gr import firdes
from optparse import OptionParser
from pss_corr import *
from sss_equ import *  
from sss_ml_fd import *
from symbol_source import *
from lte_dl_ss_source import *   
import logging

class lte_cell_scan:
	"""
	Cell scanner for LTE signals
	
	concatenate((
	  zeros((160+144*4+2048*5)/decim), gen_sss_td(N_id_1, N_id_2, True, N_re).get_data(), 
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
			top.connect(decim_lowpass, gr.file_sink(gr.sizeof_gr_complex, self.dump + "_input.cfile"))
		top.run()
		self.buffer = sink.data()
		
		logging.debug("No. samples to process: {}".format(len(self.buffer))) 
		if -1 == self.avg_frames:
			self.avg_frames = int(floor(len(self.buffer) / vec_half_frame / 2))
		logging.debug("Config: decim {}, avg_frames {}, freq_corr {}".format(self.decim, self.avg_frames, self.freq_corr)) 
		self.source = gr.vector_source_c(self.buffer)
		

	def scan(self):
		res = self.correlate_pss()
		logging.debug("Found PSS correlations " + str(res))
		for (N_id_2, frame_time, peak) in res:
			ml12 = zeros(168) 
			ml21 = zeros(168)
			ml_max = (-1e6, 0, 0)
			for N_id_1 in range(0,168):
				for slot_0_10 in range(0,1): 
					ml12[N_id_1] = self.correlate_sss(N_id_1, N_id_2, 0, frame_time, peak)
					ml21[N_id_1] = self.correlate_sss(N_id_1, N_id_2, 1, frame_time, peak)
					if (ml_max[0] < ml12[N_id_1]):
						ml_max = (ml12[N_id_1], N_id_1, 0) 
					if (ml_max[0] < ml21[N_id_1]):
						ml_max = (ml21[N_id_1], N_id_1, 1) 
			N_id_1 = ml_max[1]
			ml = ml_max[0]
			logging.debug("cell ID {}, SSS max-ml {}".format(3*N_id_1+N_id_2, ml))


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
		maxi = (0,0,0)
		for N_id_2 in range(0,3):
			pss_res = sink_pss[N_id_2].data()
			for i in range(0, len(pss_res)/2): 
				ret += [(N_id_2, pss_res[2*i], pss_res[2*i+1])]
				if maxi[2] < pss_res[2*i+1]:
					maxi = (N_id_2, pss_res[2*i], pss_res[2*i+1])
		#return ret
		return [maxi]
			

	def correlate_sss(self, N_id_1, N_id_2, slot_0_10, pss_peak_time, pss_peak):
		frame_time = pss_peak_time - (160+144*6+2048*7)
		if frame_time < 0:
			frame_time += 30720*5
		#debug_string = "Searching for SSS: N_id_1 {}, N_id_2 {}, slot_0_10 {}, frame time {:5.0f}, corr peak {:7.0f}"
		#logging.debug(debug_string.format(N_id_1, N_id_2, slot_0_10, frame_time, pss_peak))
		
		fft_size = 2048/self.decim
		N_re = 62
		top = gr.top_block("pss corr graph")
		symbol_mask = numpy.zeros(20*7)
		symbol_mask[5:7] = 1
		symbol_mask[75:77] = 1
		source = symbol_source(self.buffer, self.decim, symbol_mask, frame_time, vlen=fft_size)
		fft = gr.fft_vcc(fft_size, True, (window.blackmanharris(1024)), True, 1)
		vector_to_stream_0 = gr.vector_to_stream(gr.sizeof_gr_complex*1, fft_size)
		keep_m_in_n_0 = gr.keep_m_in_n(gr.sizeof_gr_complex, N_re, fft_size, (fft_size-N_re)/2)
		stream_to_vector_0_0 = gr.stream_to_vector(gr.sizeof_gr_complex*1, N_re)
		deinterleave_0 = gr.deinterleave(gr.sizeof_gr_complex*N_re)
		sss_equ_0 = sss_equ(N_id_2=N_id_2)
		sss_corr_0 = sss_ml_fd(decim=self.decim,N_id_1=N_id_1,N_id_2=N_id_2,avg_frames=1,slot_0_10=slot_0_10)
		ml_sink = gr.vector_sink_f(1)
		top.connect(source, fft, vector_to_stream_0, keep_m_in_n_0, stream_to_vector_0_0, deinterleave_0)
		top.connect((deinterleave_0,0), (sss_equ_0,1))
		top.connect((deinterleave_0,1), (sss_equ_0,0))
		top.connect(sss_equ_0, sss_corr_0, ml_sink)

		if self.dump != None:
			top.connect(source, gr.file_sink(gr.sizeof_gr_complex*fft_size, self.dump + "_pss{}_sss_in.cfile".format(N_id_2)))
			top.connect((deinterleave_0,0), gr.file_sink(gr.sizeof_gr_complex*N_re, self.dump + "_pss{}_sss_fd_in.cfile".format(N_id_2)))
			top.connect(sss_equ_0, gr.file_sink(gr.sizeof_gr_complex*N_re, self.dump + "_pss{}_sss_fd_equ.cfile".format(N_id_2)))
			top.connect(sss_corr_0, gr.file_sink(gr.sizeof_float, self.dump + "_pss{}_sss_corr.cfile".format(N_id_2)))

		top.run()
		return mean(ml_sink.data())


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
	#logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)
	logging.basicConfig(format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', level=logging.DEBUG)
	logging.getLogger('symbol_source').setLevel(logging.WARN)
	logging.getLogger('gen_sss_fd').setLevel(logging.WARN)
	tb = lte_cell_scan(decim=options.decim, avg_frames=options.avg_frames, freq_corr=options.freq_corr, file_name=options.file_name, dump=options.dump)
	tb.scan()


					
					