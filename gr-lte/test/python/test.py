#!/usr/bin/python

from pss_source import *
from pss_corr import *


#pss = pss_source_fd(0)
#pss = pss_source_td(0)

pss = gen_pss_td(2, 128)
#print pss.get_data()
print pss.get_data_conj_rev()

pss_corr(0, 16, 1)