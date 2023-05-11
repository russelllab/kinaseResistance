#!/usr/bin/env/Python 3.6.8
import sys
import os
import logging

script = "create_svg_20230511_kinases_GS.py"
mutinfos =  "Mutational_Infofile_Kinases_V5.txt"
konserve = "GenerelleKonservierung_May-11-2023.txt"
seqident =  "SeqIdentity_Matrix_May-11-2023.txt"


try:
	os.system("python3 "+script+" P15056 600 30 10 humanKinasesTrimmed.clustal "+mutinfos+" "+konserve+" none "+seqident+" 1")
except:
	logging.exception("message")
	pass

try:
	os.system("python3 "+script+" P15056 600 30 10 humanKinasesTrimmed.clustal "+mutinfos+" "+konserve+" none "+seqident+" 2")
except:
	logging.exception("message")
	pass

try:
	os.system("python3 "+script+" P00533 790 20 10 humanKinasesTrimmed.clustal "+mutinfos+" "+konserve+" none "+seqident+" 1")
except:
	logging.exception("message")
	pass

try:
	os.system("python3 "+script+" P00533 790 30 10 humanKinasesTrimmed.clustal "+mutinfos+" "+konserve+" none "+seqident+" 2")
except:
	logging.exception("message")
	pass

try:
	os.system("python3 "+script+" P46734 84 30 10 humanKinasesTrimmed.clustal "+mutinfos+" "+konserve+" none "+seqident+" 1")
except:
	logging.exception("message")
	pass

try:
	os.system("python3 "+script+" P46734 84 15 5 humanKinasesTrimmed.clustal "+mutinfos+" "+konserve+" none "+seqident+" 2")
except:
	logging.exception("message")
	pass

try:
	os.system("python3 "+script+" P46734 84 30000 30000 humanKinasesTrimmed.clustal "+mutinfos+" "+konserve+" none "+seqident+" 5")
except:
	logging.exception("message")
	pass


try:
	os.system("python3 "+script+" P46734 84 5 5 humanKinasesTrimmed.clustal "+mutinfos+" "+konserve+" none "+seqident+" 10")
except:
	logging.exception("message")
	pass
try:
	os.system("python3 "+script+" P46734 84 50 50 humanKinasesTrimmed.clustal "+mutinfos+" "+konserve+" none "+seqident+" 3")
except:
	logging.exception("message")
	pass
try:
	os.system("python3 "+script+" P15056 600 1 1 humanKinasesTrimmed.clustal "+mutinfos+" "+konserve+" none "+seqident+" 2")
except:
	logging.exception("message")
	pass



