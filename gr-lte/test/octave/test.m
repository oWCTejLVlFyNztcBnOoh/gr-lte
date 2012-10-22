addpath ('/home/user/projects/openlte-code/octave')
addpath("/home/user/projects/gnuradio/gnuradio-core/src/utils/")
in=read_complex_binary('/home/user/git/gr-lte/gr-lte/test/traces/lte_02_796m_30720k_frame.cfile');
% frame start: 301104 Ts
% PSS: 7th symbol 1st slot: 7072 Ts
lte_fdd_dl_receive(c)

% Generate PSS
pss0 = lte_generate_pss(0);
pss0_fd_128 = [0 pss0(1:31) zeros(1,65) pss0(32:62)];
pss0_td_128 = ifft(pss0_fd_128);
pss0_td_128 = [pss0_td_128(120:128) pss0_td_128];

pss1 = lte_generate_pss(1);
pss1_fd_128 = [0 pss1(1:31) zeros(1,65) pss1(32:62)];
pss1_td_128 = ifft(pss1_fd_128);
pss1_td_128 = [pss1_td_128(120:128) pss1_td_128];

pss2 = lte_generate_pss(2);
pss2_fd_128 = [0 pss2(1:31) zeros(1,65) pss2(32:62)];
pss2_td_128 = ifft(pss2_fd_128);
pss2_td_128 = [pss2_td_128(120:128) pss2_td_128];

%pss_corr = reshape(abs(filter((pss0_td_128),1,in(1:16:30720*100))),307200/16,10);
%pss_corr = reshape(abs(filter(conj(pss0_td_128),1,in(1:16:30720*100))),307200/16,10);
pss_corr = reshape(abs(filter(conj(rot90(pss0_td_128,2)),1,in(1:16:30720*100))),307200/16,10);
plot([1:16:307200]-(2048+144)*0,sum(filter(ones(1,3)/3,1,pss_corr),2))



start=1;
endl=30720;
offset=2048;
plot(abs(filter(ones(144,1)/144,1, conj(in(start:endl)).*in(offset+start:endl+offset))));
