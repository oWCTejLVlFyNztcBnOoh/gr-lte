function [pss_corr, pss_corr_add, in] = pss_plot(prefix='foo', data_path='', decim=16)
addpath ('/home/user/projects/openlte-code/octave')
addpath("/home/user/projects/gnuradio/gnuradio-core/src/utils/")

%if nargin < 1
%    prefix = 'foo';
%end
%if nargin < 2
%    data_path = '';
%end 

in = read_complex_binary([data_path prefix '_input.cfile']);
pss_corr = read_float_binary([data_path prefix '_pss_corr_f.cfile']);
pss_corr_add = read_float_binary([data_path prefix '_pss_corr_add_f.cfile']);

idx = [0:length(pss_corr_add)-1]*decim;
plot(idx, pss_corr_add);

figure;plot(idx,pss_corr_add/max(pss_corr_add), idx,abs(in(1:length(pss_corr_add))));

sss0_fd = read_complex_binary([data_path prefix '_pss0_sss_fd.cfile']);

figure; plot(sss0_fd(1:128), 'x');
figure; plot(sss0_fd(129:256), 'x');

np/(1+npi/np)
1/(np+npi)