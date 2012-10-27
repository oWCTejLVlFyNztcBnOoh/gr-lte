#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Sss Corr3 Gui
# Generated: Sat Oct 27 19:36:17 2012
##################################################

execfile("/home/user/.grc_gnuradio/pss_chan_est2.py")
from gnuradio import digital
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio import window
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.gr import firdes
from gnuradio.wxgui import forms
from gnuradio.wxgui import scopesink2
from grc_gnuradio import blks2 as grc_blks2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
from pss_corr import *  
from sss_corr import *   
import wx

class sss_corr3_gui(grc_wxgui.top_block_gui):

	def __init__(self, freq_corr=0, avg_frames=1, N_id_2=0, N_id_1=134, decim=16):
		grc_wxgui.top_block_gui.__init__(self, title="Sss Corr3 Gui")
		_icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
		self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

		##################################################
		# Parameters
		##################################################
		self.freq_corr = freq_corr
		self.avg_frames = avg_frames
		self.N_id_2 = N_id_2
		self.N_id_1 = N_id_1
		self.decim = decim

		##################################################
		# Variables
		##################################################
		self.vec_half_frame = vec_half_frame = 30720*5/decim
		self.symbol_start = symbol_start = 144/decim
		self.samp_rate = samp_rate = 30720e3/decim
		self.rot = rot = 0
		self.noise_level = noise_level = 0
		self.fft_size = fft_size = 2048/decim
		self.N_re = N_re = 62

		##################################################
		# Blocks
		##################################################
		_symbol_start_sizer = wx.BoxSizer(wx.VERTICAL)
		self._symbol_start_text_box = forms.text_box(
			parent=self.GetWin(),
			sizer=_symbol_start_sizer,
			value=self.symbol_start,
			callback=self.set_symbol_start,
			label='symbol_start',
			converter=forms.int_converter(),
			proportion=0,
		)
		self._symbol_start_slider = forms.slider(
			parent=self.GetWin(),
			sizer=_symbol_start_sizer,
			value=self.symbol_start,
			callback=self.set_symbol_start,
			minimum=0,
			maximum=144/decim,
			num_steps=144/decim,
			style=wx.SL_HORIZONTAL,
			cast=int,
			proportion=1,
		)
		self.Add(_symbol_start_sizer)
		_rot_sizer = wx.BoxSizer(wx.VERTICAL)
		self._rot_text_box = forms.text_box(
			parent=self.GetWin(),
			sizer=_rot_sizer,
			value=self.rot,
			callback=self.set_rot,
			label='rot',
			converter=forms.float_converter(),
			proportion=0,
		)
		self._rot_slider = forms.slider(
			parent=self.GetWin(),
			sizer=_rot_sizer,
			value=self.rot,
			callback=self.set_rot,
			minimum=0,
			maximum=1,
			num_steps=100,
			style=wx.SL_HORIZONTAL,
			cast=float,
			proportion=1,
		)
		self.Add(_rot_sizer)
		self.notebook_0 = self.notebook_0 = wx.Notebook(self.GetWin(), style=wx.NB_TOP)
		self.notebook_0.AddPage(grc_wxgui.Panel(self.notebook_0), "SSS ML")
		self.notebook_0.AddPage(grc_wxgui.Panel(self.notebook_0), "SSS fd est/ref")
		self.notebook_0.AddPage(grc_wxgui.Panel(self.notebook_0), "SSS est-ref diff")
		self.notebook_0.AddPage(grc_wxgui.Panel(self.notebook_0), "input td")
		self.Add(self.notebook_0)
		_noise_level_sizer = wx.BoxSizer(wx.VERTICAL)
		self._noise_level_text_box = forms.text_box(
			parent=self.GetWin(),
			sizer=_noise_level_sizer,
			value=self.noise_level,
			callback=self.set_noise_level,
			label='noise_level',
			converter=forms.float_converter(),
			proportion=0,
		)
		self._noise_level_slider = forms.slider(
			parent=self.GetWin(),
			sizer=_noise_level_sizer,
			value=self.noise_level,
			callback=self.set_noise_level,
			minimum=0,
			maximum=10,
			num_steps=100,
			style=wx.SL_HORIZONTAL,
			cast=float,
			proportion=1,
		)
		self.Add(_noise_level_sizer)
		self.wxgui_scopesink2_0_1_1 = scopesink2.scope_sink_f(
			self.notebook_0.GetPage(0).GetWin(),
			title="Scope Plot",
			sample_rate=samp_rate/N_re,
			v_scale=0,
			v_offset=0,
			t_scale=0,
			ac_couple=False,
			xy_mode=False,
			num_inputs=1,
			trig_mode=gr.gr_TRIG_MODE_AUTO,
			y_axis_label="Counts",
		)
		self.notebook_0.GetPage(0).Add(self.wxgui_scopesink2_0_1_1.win)
		self.wxgui_scopesink2_0_1_0_0 = scopesink2.scope_sink_c(
			self.notebook_0.GetPage(3).GetWin(),
			title="Scope Plot",
			sample_rate=samp_rate,
			v_scale=0,
			v_offset=0,
			t_scale=0,
			ac_couple=False,
			xy_mode=False,
			num_inputs=1,
			trig_mode=gr.gr_TRIG_MODE_AUTO,
			y_axis_label="Counts",
		)
		self.notebook_0.GetPage(3).Add(self.wxgui_scopesink2_0_1_0_0.win)
		self.wxgui_scopesink2_0_1_0 = scopesink2.scope_sink_c(
			self.notebook_0.GetPage(1).GetWin(),
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
		self.notebook_0.GetPage(1).Add(self.wxgui_scopesink2_0_1_0.win)
		self.wxgui_scopesink2_0_1 = scopesink2.scope_sink_f(
			self.notebook_0.GetPage(2).GetWin(),
			title="Scope Plot",
			sample_rate=samp_rate,
			v_scale=0,
			v_offset=0,
			t_scale=0,
			ac_couple=False,
			xy_mode=False,
			num_inputs=3,
			trig_mode=gr.gr_TRIG_MODE_AUTO,
			y_axis_label="Counts",
		)
		self.notebook_0.GetPage(2).Add(self.wxgui_scopesink2_0_1.win)
		self.pss_chan_est2_0 = pss_chan_est2(
			N_id_2=N_id_2,
		)
		self.gr_vector_to_stream_1 = gr.vector_to_stream(gr.sizeof_gr_complex*1, fft_size)
		self.gr_vector_to_stream_0_2_0_1 = gr.vector_to_stream(gr.sizeof_gr_complex*1, N_re)
		self.gr_vector_to_stream_0_2_0_0_0 = gr.vector_to_stream(gr.sizeof_gr_complex*1, N_re)
		self.gr_vector_to_stream_0_2_0_0 = gr.vector_to_stream(gr.sizeof_gr_complex*1, N_re)
		self.gr_vector_to_stream_0_2_0 = gr.vector_to_stream(gr.sizeof_gr_complex*1, N_re)
		self.gr_vector_to_stream_0_2 = gr.vector_to_stream(gr.sizeof_float*1, N_re)
		self.gr_vector_to_stream_0_0 = gr.vector_to_stream(gr.sizeof_gr_complex*1, N_re)
		self.gr_vector_to_stream_0 = gr.vector_to_stream(gr.sizeof_gr_complex*1, fft_size)
		self.gr_vector_source_x_0_1 = gr.vector_source_c((gen_sss_fd( N_id_1, N_id_2, fft_size).get_sss(False)), True, fft_size)
		self.gr_vector_source_x_0_0_1 = gr.vector_source_c((gen_pss_fd(N_id_2, fft_size, False).get_data()), True, fft_size)
		self.gr_vector_source_x_0_0_0_0 = gr.vector_source_c((numpy.conjugate(gen_sss_fd(N_id_1,N_id_2, N_re).get_sss(False))), True, N_re)
		self.gr_vector_source_x_0_0_0 = gr.vector_source_c((numpy.conjugate(gen_sss_fd(N_id_1,N_id_2, N_re).get_sss(True))), True, N_re)
		self.gr_vector_source_x_0_0 = gr.vector_source_c((gen_pss_fd(N_id_2, fft_size, False).get_data()), True, fft_size)
		self.gr_vector_source_x_0 = gr.vector_source_c((gen_sss_fd( N_id_1, N_id_2, fft_size).get_sss(True)), True, fft_size)
		self.gr_throttle_0 = gr.throttle(gr.sizeof_gr_complex*1, samp_rate)
		self.gr_sub_xx_0 = gr.sub_cc(N_re)
		self.gr_stream_to_vector_1_0 = gr.stream_to_vector(gr.sizeof_gr_complex*1, N_re)
		self.gr_stream_to_vector_0_0 = gr.stream_to_vector(gr.sizeof_gr_complex*1, N_re)
		self.gr_stream_to_vector_0 = gr.stream_to_vector(gr.sizeof_gr_complex*1, fft_size)
		self.gr_repeat_0 = gr.repeat(gr.sizeof_gr_complex*1, N_re)
		self.gr_noise_source_x_0 = gr.noise_source_c(gr.GR_GAUSSIAN, noise_level, 0)
		self.gr_multiply_xx_2 = gr.multiply_vcc(N_re)
		self.gr_multiply_xx_1 = gr.multiply_vcc(N_re)
		self.gr_multiply_xx_0 = gr.multiply_vcc(N_re)
		self.gr_multiply_const_vxx_1 = gr.multiply_const_vff((-1e-0, ))
		self.gr_multiply_const_vxx_0 = gr.multiply_const_vcc((0.005*exp(rot*2*numpy.pi*1j), ))
		self.gr_keep_m_in_n_1 = gr.keep_m_in_n(gr.sizeof_gr_complex, fft_size, fft_size+144/decim, symbol_start)
		self.gr_keep_m_in_n_0 = gr.keep_m_in_n(gr.sizeof_gr_complex, N_re, fft_size, (fft_size-N_re)/2)
		self.gr_interleave_1 = gr.interleave(gr.sizeof_gr_complex*N_re)
		self.gr_interleave_0 = gr.interleave(gr.sizeof_gr_complex*fft_size)
		self.gr_file_source_0 = gr.file_source(gr.sizeof_gr_complex*1, "/home/user/git/gr-lte/gr-lte/test/foo_pss0_sss_in.cfile", True)
		self.gr_fft_vxx_1 = gr.fft_vcc(fft_size, False, (window.blackmanharris(1024)), True, 1)
		self.gr_fft_vxx_0 = gr.fft_vcc(fft_size, True, (window.blackmanharris(1024)), True, 1)
		self.gr_deinterleave_0 = gr.deinterleave(gr.sizeof_gr_complex*N_re)
		self.gr_conjugate_cc_1 = gr.conjugate_cc()
		self.gr_complex_to_mag_squared_0 = gr.complex_to_mag_squared(N_re)
		self.gr_complex_to_float_0 = gr.complex_to_float(1)
		self.gr_channel_model_0 = gr.channel_model(
			noise_voltage=0.005*noise_level,
			frequency_offset=0.0,
			epsilon=1,
			taps=(0.005*exp(rot*2*numpy.pi*1j), ),
			noise_seed=0,
		)
		self.gr_add_xx_0 = gr.add_vcc(1)
		self.fir_filter_xxx_0_0 = filter.fir_filter_fff(N_re, (numpy.ones(N_re)/N_re))
		self.fir_filter_xxx_0 = filter.fir_filter_ccc(N_re, (numpy.ones(N_re)/N_re))
		self.digital_ofdm_cyclic_prefixer_0 = digital.ofdm_cyclic_prefixer(fft_size, fft_size+144/decim)
		self.blks2_selector_0_0 = grc_blks2.selector(
			item_size=gr.sizeof_gr_complex*1,
			num_inputs=2,
			num_outputs=1,
			input_index=1,
			output_index=0,
		)
		self.blks2_selector_0 = grc_blks2.selector(
			item_size=gr.sizeof_gr_complex*1,
			num_inputs=3,
			num_outputs=2,
			input_index=0,
			output_index=1,
		)

		##################################################
		# Connections
		##################################################
		self.connect((self.gr_stream_to_vector_0, 0), (self.gr_fft_vxx_0, 0))
		self.connect((self.gr_deinterleave_0, 1), (self.pss_chan_est2_0, 0))
		self.connect((self.gr_stream_to_vector_0_0, 0), (self.gr_deinterleave_0, 0))
		self.connect((self.gr_keep_m_in_n_0, 0), (self.gr_stream_to_vector_0_0, 0))
		self.connect((self.gr_vector_to_stream_0, 0), (self.gr_keep_m_in_n_0, 0))
		self.connect((self.gr_fft_vxx_0, 0), (self.gr_vector_to_stream_0, 0))
		self.connect((self.pss_chan_est2_0, 0), (self.gr_multiply_xx_1, 1))
		self.connect((self.gr_deinterleave_0, 0), (self.gr_multiply_xx_1, 0))
		self.connect((self.gr_vector_source_x_0_0_0, 0), (self.gr_interleave_1, 0))
		self.connect((self.gr_vector_source_x_0_0_0_0, 0), (self.gr_interleave_1, 1))
		self.connect((self.gr_interleave_1, 0), (self.gr_multiply_xx_0, 0))
		self.connect((self.gr_multiply_xx_0, 0), (self.gr_vector_to_stream_0_0, 0))
		self.connect((self.gr_multiply_xx_2, 0), (self.gr_sub_xx_0, 0))
		self.connect((self.gr_stream_to_vector_1_0, 0), (self.gr_multiply_xx_2, 1))
		self.connect((self.gr_interleave_1, 0), (self.gr_multiply_xx_2, 0))
		self.connect((self.gr_multiply_xx_1, 0), (self.gr_sub_xx_0, 1))
		self.connect((self.gr_sub_xx_0, 0), (self.gr_complex_to_mag_squared_0, 0))
		self.connect((self.gr_complex_to_mag_squared_0, 0), (self.gr_vector_to_stream_0_2, 0))
		self.connect((self.gr_vector_to_stream_0_2, 0), (self.gr_multiply_const_vxx_1, 0))
		self.connect((self.gr_multiply_xx_1, 0), (self.gr_multiply_xx_0, 1))
		self.connect((self.gr_vector_to_stream_0_0, 0), (self.fir_filter_xxx_0, 0))
		self.connect((self.fir_filter_xxx_0, 0), (self.gr_conjugate_cc_1, 0))
		self.connect((self.gr_conjugate_cc_1, 0), (self.gr_repeat_0, 0))
		self.connect((self.gr_repeat_0, 0), (self.gr_stream_to_vector_1_0, 0))
		self.connect((self.gr_multiply_xx_1, 0), (self.gr_vector_to_stream_0_2_0, 0))
		self.connect((self.gr_multiply_const_vxx_1, 0), (self.wxgui_scopesink2_0_1, 2))
		self.connect((self.gr_complex_to_float_0, 1), (self.wxgui_scopesink2_0_1, 1))
		self.connect((self.gr_complex_to_float_0, 0), (self.wxgui_scopesink2_0_1, 0))
		self.connect((self.gr_vector_to_stream_0_2_0, 0), (self.gr_complex_to_float_0, 0))
		self.connect((self.fir_filter_xxx_0_0, 0), (self.wxgui_scopesink2_0_1_1, 0))
		self.connect((self.gr_multiply_const_vxx_1, 0), (self.fir_filter_xxx_0_0, 0))
		self.connect((self.gr_interleave_1, 0), (self.gr_vector_to_stream_0_2_0_0_0, 0))
		self.connect((self.gr_vector_to_stream_0_2_0_0_0, 0), (self.wxgui_scopesink2_0_1_0, 0))
		self.connect((self.gr_multiply_xx_1, 0), (self.gr_vector_to_stream_0_2_0_0, 0))
		self.connect((self.gr_vector_to_stream_0_2_0_0, 0), (self.wxgui_scopesink2_0_1_0, 1))
		self.connect((self.gr_noise_source_x_0, 0), (self.gr_add_xx_0, 1))
		self.connect((self.gr_add_xx_0, 0), (self.gr_multiply_const_vxx_0, 0))
		self.connect((self.gr_throttle_0, 0), (self.gr_stream_to_vector_0, 0))
		self.connect((self.blks2_selector_0, 0), (self.gr_channel_model_0, 0))
		self.connect((self.blks2_selector_0, 1), (self.gr_add_xx_0, 0))
		self.connect((self.gr_channel_model_0, 0), (self.blks2_selector_0_0, 0))
		self.connect((self.gr_multiply_const_vxx_0, 0), (self.blks2_selector_0_0, 1))
		self.connect((self.gr_keep_m_in_n_1, 0), (self.blks2_selector_0, 2))
		self.connect((self.gr_vector_to_stream_1, 0), (self.blks2_selector_0, 1))
		self.connect((self.gr_fft_vxx_1, 0), (self.digital_ofdm_cyclic_prefixer_0, 0))
		self.connect((self.digital_ofdm_cyclic_prefixer_0, 0), (self.gr_keep_m_in_n_1, 0))
		self.connect((self.gr_file_source_0, 0), (self.blks2_selector_0, 0))
		self.connect((self.gr_fft_vxx_1, 0), (self.gr_vector_to_stream_1, 0))
		self.connect((self.gr_interleave_0, 0), (self.gr_fft_vxx_1, 0))
		self.connect((self.blks2_selector_0_0, 0), (self.gr_throttle_0, 0))
		self.connect((self.gr_vector_source_x_0_0, 0), (self.gr_interleave_0, 1))
		self.connect((self.gr_vector_source_x_0_0_1, 0), (self.gr_interleave_0, 3))
		self.connect((self.gr_vector_source_x_0, 0), (self.gr_interleave_0, 0))
		self.connect((self.gr_vector_source_x_0_1, 0), (self.gr_interleave_0, 2))
		self.connect((self.gr_vector_to_stream_0_2_0_1, 0), (self.wxgui_scopesink2_0_1_0_0, 0))
		self.connect((self.gr_deinterleave_0, 1), (self.gr_vector_to_stream_0_2_0_1, 0))

	def get_freq_corr(self):
		return self.freq_corr

	def set_freq_corr(self, freq_corr):
		self.freq_corr = freq_corr

	def get_avg_frames(self):
		return self.avg_frames

	def set_avg_frames(self, avg_frames):
		self.avg_frames = avg_frames

	def get_N_id_2(self):
		return self.N_id_2

	def set_N_id_2(self, N_id_2):
		self.N_id_2 = N_id_2
		self.pss_chan_est2_0.set_N_id_2(self.N_id_2)

	def get_N_id_1(self):
		return self.N_id_1

	def set_N_id_1(self, N_id_1):
		self.N_id_1 = N_id_1

	def get_decim(self):
		return self.decim

	def set_decim(self, decim):
		self.decim = decim
		self.set_vec_half_frame(30720*5/self.decim)
		self.set_symbol_start(144/self.decim)
		self.set_samp_rate(30720e3/self.decim)
		self.set_fft_size(2048/self.decim)
		self.gr_keep_m_in_n_1.set_n(self.fft_size+144/self.decim)

	def get_vec_half_frame(self):
		return self.vec_half_frame

	def set_vec_half_frame(self, vec_half_frame):
		self.vec_half_frame = vec_half_frame

	def get_symbol_start(self):
		return self.symbol_start

	def set_symbol_start(self, symbol_start):
		self.symbol_start = symbol_start
		self._symbol_start_slider.set_value(self.symbol_start)
		self._symbol_start_text_box.set_value(self.symbol_start)
		self.gr_keep_m_in_n_1.set_offset(self.symbol_start)

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate
		self.wxgui_scopesink2_0_1_1.set_sample_rate(self.samp_rate/self.N_re)
		self.wxgui_scopesink2_0_1.set_sample_rate(self.samp_rate)
		self.wxgui_scopesink2_0_1_0.set_sample_rate(self.samp_rate)
		self.wxgui_scopesink2_0_1_0_0.set_sample_rate(self.samp_rate)

	def get_rot(self):
		return self.rot

	def set_rot(self, rot):
		self.rot = rot
		self._rot_slider.set_value(self.rot)
		self._rot_text_box.set_value(self.rot)
		self.gr_channel_model_0.set_taps((0.005*exp(self.rot*2*numpy.pi*1j), ))
		self.gr_multiply_const_vxx_0.set_k((0.005*exp(self.rot*2*numpy.pi*1j), ))

	def get_noise_level(self):
		return self.noise_level

	def set_noise_level(self, noise_level):
		self.noise_level = noise_level
		self.gr_channel_model_0.set_noise_voltage(0.005*self.noise_level)
		self._noise_level_slider.set_value(self.noise_level)
		self._noise_level_text_box.set_value(self.noise_level)
		self.gr_noise_source_x_0.set_amplitude(self.noise_level)

	def get_fft_size(self):
		return self.fft_size

	def set_fft_size(self, fft_size):
		self.fft_size = fft_size
		self.gr_keep_m_in_n_0.set_offset((self.fft_size-self.N_re)/2)
		self.gr_keep_m_in_n_0.set_n(self.fft_size)
		self.gr_keep_m_in_n_1.set_m(self.fft_size)
		self.gr_keep_m_in_n_1.set_n(self.fft_size+144/self.decim)

	def get_N_re(self):
		return self.N_re

	def set_N_re(self, N_re):
		self.N_re = N_re
		self.gr_keep_m_in_n_0.set_offset((self.fft_size-self.N_re)/2)
		self.gr_keep_m_in_n_0.set_m(self.N_re)
		self.fir_filter_xxx_0.set_taps((numpy.ones(self.N_re)/self.N_re))
		self.fir_filter_xxx_0_0.set_taps((numpy.ones(self.N_re)/self.N_re))
		self.wxgui_scopesink2_0_1_1.set_sample_rate(self.samp_rate/self.N_re)

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	parser.add_option("", "--freq-corr", dest="freq_corr", type="eng_float", default=eng_notation.num_to_str(0),
		help="Set freq_corr [default=%default]")
	parser.add_option("", "--avg-frames", dest="avg_frames", type="intx", default=1,
		help="Set avg_frames [default=%default]")
	parser.add_option("", "--N-id-2", dest="N_id_2", type="intx", default=0,
		help="Set N_id_2 [default=%default]")
	parser.add_option("", "--N-id-1", dest="N_id_1", type="intx", default=134,
		help="Set N_id_1 [default=%default]")
	parser.add_option("", "--decim", dest="decim", type="intx", default=16,
		help="Set decim [default=%default]")
	(options, args) = parser.parse_args()
	tb = sss_corr3_gui(freq_corr=options.freq_corr, avg_frames=options.avg_frames, N_id_2=options.N_id_2, N_id_1=options.N_id_1, decim=options.decim)
	tb.Run(True)

