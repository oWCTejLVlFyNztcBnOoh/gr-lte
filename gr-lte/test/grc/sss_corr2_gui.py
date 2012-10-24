#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Sss Corr2 Gui
# Generated: Wed Oct 24 23:29:20 2012
##################################################

from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import window
from gnuradio.eng_option import eng_option
from gnuradio.gr import firdes
from gnuradio.wxgui import forms
from gnuradio.wxgui import scopesink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
from pss_corr import *  
from sss_corr import *   
import wx

class sss_corr2_gui(grc_wxgui.top_block_gui):

	def __init__(self, freq_corr=0, N_id_1=134, avg_frames=1, N_id_2=0, decim=16):
		grc_wxgui.top_block_gui.__init__(self, title="Sss Corr2 Gui")
		_icon_path = "/usr/share/icons/hicolor/32x32/apps/gnuradio-grc.png"
		self.SetIcon(wx.Icon(_icon_path, wx.BITMAP_TYPE_ANY))

		##################################################
		# Parameters
		##################################################
		self.freq_corr = freq_corr
		self.N_id_1 = N_id_1
		self.avg_frames = avg_frames
		self.N_id_2 = N_id_2
		self.decim = decim

		##################################################
		# Variables
		##################################################
		self.vec_half_frame = vec_half_frame = 30720*5/decim
		self.samp_rate = samp_rate = 30720e3/decim
		self.rot = rot = 0
		self.noise_level = noise_level = 0
		self.fft_size = fft_size = 2048/decim

		##################################################
		# Blocks
		##################################################
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
		self.wxgui_scopesink2_0 = scopesink2.scope_sink_c(
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
		self.gr_vector_to_stream_1 = gr.vector_to_stream(gr.sizeof_gr_complex*1, fft_size)
		self.gr_vector_to_stream_0_2 = gr.vector_to_stream(gr.sizeof_gr_complex*1, fft_size)
		self.gr_vector_to_stream_0_1 = gr.vector_to_stream(gr.sizeof_gr_complex*1, fft_size)
		self.gr_vector_to_stream_0_0 = gr.vector_to_stream(gr.sizeof_gr_complex*1, fft_size)
		self.gr_vector_to_stream_0 = gr.vector_to_stream(gr.sizeof_float*1, fft_size)
		self.gr_vector_source_x_0_0_0 = gr.vector_source_c((gen_pss_fd(N_id_2, fft_size, False).get_data()), True, fft_size)
		self.gr_vector_source_x_0_0 = gr.vector_source_c((gen_pss_fd(N_id_2, fft_size, False).get_data()), True, fft_size)
		self.gr_vector_source_x_0 = gr.vector_source_c((gen_sss_fd( N_id_1, N_id_2, fft_size).get_sss(True)), True, fft_size)
		self.gr_stream_to_vector_0_0 = gr.stream_to_vector(gr.sizeof_gr_complex*1, fft_size)
		self.gr_stream_to_vector_0 = gr.stream_to_vector(gr.sizeof_gr_complex*1, fft_size)
		self.gr_repeat_0 = gr.repeat(gr.sizeof_float*1, fft_size)
		self.gr_null_sink_0 = gr.null_sink(gr.sizeof_gr_complex*1)
		self.gr_noise_source_x_0 = gr.noise_source_c(gr.GR_GAUSSIAN, noise_level, 0)
		self.gr_multiply_xx_1 = gr.multiply_vcc(1)
		self.gr_multiply_xx_0 = gr.multiply_vcc(fft_size)
		self.gr_multiply_const_vxx_1 = gr.multiply_const_vcc((1/1500., ))
		self.gr_multiply_const_vxx_0 = gr.multiply_const_vcc((exp(rot*2*numpy.pi*1j), ))
		self.gr_interleave_0 = gr.interleave(gr.sizeof_gr_complex*fft_size)
		self.gr_integrate_xx_0 = gr.integrate_ff(fft_size)
		self.gr_float_to_complex_0_0 = gr.float_to_complex(1)
		self.gr_float_to_complex_0 = gr.float_to_complex(1)
		self.gr_fft_vxx_1 = gr.fft_vcc(fft_size, False, (window.blackmanharris(1024)), True, 1)
		self.gr_fft_vxx_0 = gr.fft_vcc(fft_size, True, (window.blackmanharris(1024)), True, 1)
		self.gr_divide_xx_0_1 = gr.divide_cc(1)
		self.gr_divide_xx_0_0 = gr.divide_ff(1)
		self.gr_divide_xx_0 = gr.divide_cc(1)
		self.gr_deinterleave_0 = gr.deinterleave(gr.sizeof_gr_complex*fft_size)
		self.gr_conjugate_cc_1 = gr.conjugate_cc()
		self.gr_conjugate_cc_0 = gr.conjugate_cc()
		self.gr_complex_to_mag_squared_0_0 = gr.complex_to_mag_squared(1)
		self.gr_complex_to_mag_squared_0 = gr.complex_to_mag_squared(fft_size)
		self.gr_add_xx_0 = gr.add_vcc(1)
		self.gr_add_const_vxx_0 = gr.add_const_vff((1, ))
		self.const_source_x_0_0 = gr.sig_source_f(0, gr.GR_CONST_WAVE, 0, 0, 0)
		self.const_source_x_0 = gr.sig_source_f(0, gr.GR_CONST_WAVE, 0, 0, 0)

		##################################################
		# Connections
		##################################################
		self.connect((self.gr_conjugate_cc_0, 0), (self.gr_stream_to_vector_0_0, 0))
		self.connect((self.gr_stream_to_vector_0_0, 0), (self.gr_multiply_xx_0, 1))
		self.connect((self.gr_deinterleave_0, 0), (self.gr_multiply_xx_0, 0))
		self.connect((self.gr_multiply_xx_0, 0), (self.gr_complex_to_mag_squared_0, 0))
		self.connect((self.gr_complex_to_mag_squared_0, 0), (self.gr_vector_to_stream_0, 0))
		self.connect((self.gr_vector_to_stream_0, 0), (self.gr_integrate_xx_0, 0))
		self.connect((self.gr_integrate_xx_0, 0), (self.gr_repeat_0, 0))
		self.connect((self.gr_multiply_xx_0, 0), (self.gr_vector_to_stream_0_0, 0))
		self.connect((self.gr_vector_to_stream_0_1, 0), (self.gr_multiply_xx_1, 1))
		self.connect((self.gr_divide_xx_0, 0), (self.gr_multiply_xx_1, 0))
		self.connect((self.gr_float_to_complex_0, 0), (self.gr_divide_xx_0, 1))
		self.connect((self.gr_conjugate_cc_1, 0), (self.gr_divide_xx_0, 0))
		self.connect((self.gr_deinterleave_0, 1), (self.gr_vector_to_stream_0_1, 0))
		self.connect((self.gr_repeat_0, 0), (self.gr_float_to_complex_0, 0))
		self.connect((self.const_source_x_0, 0), (self.gr_float_to_complex_0, 1))
		self.connect((self.gr_vector_to_stream_0_0, 0), (self.gr_conjugate_cc_1, 0))
		self.connect((self.gr_vector_to_stream_0_0, 0), (self.gr_complex_to_mag_squared_0_0, 0))
		self.connect((self.gr_complex_to_mag_squared_0_0, 0), (self.gr_divide_xx_0_0, 0))
		self.connect((self.gr_repeat_0, 0), (self.gr_divide_xx_0_0, 1))
		self.connect((self.gr_divide_xx_0_0, 0), (self.gr_add_const_vxx_0, 0))
		self.connect((self.gr_add_const_vxx_0, 0), (self.gr_float_to_complex_0_0, 0))
		self.connect((self.const_source_x_0_0, 0), (self.gr_float_to_complex_0_0, 1))
		self.connect((self.gr_float_to_complex_0_0, 0), (self.gr_divide_xx_0_1, 1))
		self.connect((self.gr_multiply_xx_1, 0), (self.gr_divide_xx_0_1, 0))
		self.connect((self.gr_divide_xx_0_1, 0), (self.wxgui_scopesink2_0, 0))
		self.connect((self.gr_stream_to_vector_0, 0), (self.gr_fft_vxx_0, 0))
		self.connect((self.gr_fft_vxx_0, 0), (self.gr_deinterleave_0, 0))
		self.connect((self.gr_vector_to_stream_0_2, 0), (self.gr_conjugate_cc_0, 0))
		self.connect((self.gr_divide_xx_0_1, 0), (self.gr_null_sink_0, 0))
		self.connect((self.gr_vector_source_x_0, 0), (self.gr_interleave_0, 1))
		self.connect((self.gr_interleave_0, 0), (self.gr_fft_vxx_1, 0))
		self.connect((self.gr_fft_vxx_1, 0), (self.gr_vector_to_stream_1, 0))
		self.connect((self.gr_multiply_const_vxx_0, 0), (self.gr_stream_to_vector_0, 0))
		self.connect((self.gr_vector_source_x_0_0, 0), (self.gr_interleave_0, 0))
		self.connect((self.gr_vector_source_x_0_0_0, 0), (self.gr_vector_to_stream_0_2, 0))
		self.connect((self.gr_noise_source_x_0, 0), (self.gr_add_xx_0, 1))
		self.connect((self.gr_vector_to_stream_1, 0), (self.gr_add_xx_0, 0))
		self.connect((self.gr_add_xx_0, 0), (self.gr_multiply_const_vxx_0, 0))
		self.connect((self.gr_multiply_const_vxx_1, 0), (self.wxgui_scopesink2_0, 1))
		self.connect((self.gr_vector_to_stream_0_1, 0), (self.gr_multiply_const_vxx_1, 0))

	def get_freq_corr(self):
		return self.freq_corr

	def set_freq_corr(self, freq_corr):
		self.freq_corr = freq_corr

	def get_N_id_1(self):
		return self.N_id_1

	def set_N_id_1(self, N_id_1):
		self.N_id_1 = N_id_1

	def get_avg_frames(self):
		return self.avg_frames

	def set_avg_frames(self, avg_frames):
		self.avg_frames = avg_frames

	def get_N_id_2(self):
		return self.N_id_2

	def set_N_id_2(self, N_id_2):
		self.N_id_2 = N_id_2

	def get_decim(self):
		return self.decim

	def set_decim(self, decim):
		self.decim = decim
		self.set_vec_half_frame(30720*5/self.decim)
		self.set_samp_rate(30720e3/self.decim)
		self.set_fft_size(2048/self.decim)

	def get_vec_half_frame(self):
		return self.vec_half_frame

	def set_vec_half_frame(self, vec_half_frame):
		self.vec_half_frame = vec_half_frame

	def get_samp_rate(self):
		return self.samp_rate

	def set_samp_rate(self, samp_rate):
		self.samp_rate = samp_rate
		self.wxgui_scopesink2_0.set_sample_rate(self.samp_rate)

	def get_rot(self):
		return self.rot

	def set_rot(self, rot):
		self.rot = rot
		self._rot_slider.set_value(self.rot)
		self._rot_text_box.set_value(self.rot)
		self.gr_multiply_const_vxx_0.set_k((exp(self.rot*2*numpy.pi*1j), ))

	def get_noise_level(self):
		return self.noise_level

	def set_noise_level(self, noise_level):
		self.noise_level = noise_level
		self.gr_noise_source_x_0.set_amplitude(self.noise_level)
		self._noise_level_slider.set_value(self.noise_level)
		self._noise_level_text_box.set_value(self.noise_level)

	def get_fft_size(self):
		return self.fft_size

	def set_fft_size(self, fft_size):
		self.fft_size = fft_size

if __name__ == '__main__':
	parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
	parser.add_option("", "--freq-corr", dest="freq_corr", type="eng_float", default=eng_notation.num_to_str(0),
		help="Set freq_corr [default=%default]")
	parser.add_option("", "--N-id-1", dest="N_id_1", type="intx", default=134,
		help="Set N_id_1 [default=%default]")
	parser.add_option("", "--avg-frames", dest="avg_frames", type="intx", default=1,
		help="Set avg_frames [default=%default]")
	parser.add_option("", "--N-id-2", dest="N_id_2", type="intx", default=0,
		help="Set N_id_2 [default=%default]")
	parser.add_option("", "--decim", dest="decim", type="intx", default=16,
		help="Set decim [default=%default]")
	(options, args) = parser.parse_args()
	tb = sss_corr2_gui(freq_corr=options.freq_corr, N_id_1=options.N_id_1, avg_frames=options.avg_frames, N_id_2=options.N_id_2, decim=options.decim)
	tb.Run(True)

