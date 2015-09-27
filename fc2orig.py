#!/usr/local/bin/python3

# the weather data is in fc_FINAL_campanas.meteo and looks like
#
# #JulianDate	P         H 	T	         D	           W
# 2453602.307	29.967	  13  57.799999	 15	          40
#
# the seeing  is in fc_final_dimm4.seeing and looks like:
#
# #JulianDate	Airmass	na   	na	  na    	na	 na	  seeing
# 2453644.493	1.13	  11730	0.807	0.758	29091	41.4	0.737
#
# great…here are the conversions, for P, T, W.
#
# symbol	   P	     H	    T	        D	     W
# old unit	InHg	    %	    deg F	    deg	   mph
# new unit	hPa	      %	    deg C	    deg	   m/s
#          	 24.76520098	 	    (T-32)*5/9	 	   0.44704

#parameters = ['S', 'T', 'H', 'P', 'W', 'E', 'R', 'TETA0']

import glob
import os
import ambient
import utils

parser = ambient.ArgumentParser (description="""
this script only takes --indir and --outdir arguments.
it reads in --indir exactly 1 file ending in .meteo and 1
file ending in .seeing and produces in --outdir the files S.orig, T.orig, H.orig, P.orig, and W.orig
and empty files called E.orig, R.orig and TETAO.orig.
""")

#parse user's command line arguments
args = parser.parse_args ()

#make sure the output directory exists
utils.makeDir (args.outdir)

seeingPaths = glob.glob (os.path.join (args.indir,'*.seeing'))
meteoPaths = glob.glob (os.path.join (args.indir,'*.meteo'))
assert len (seeingPaths) == 1
assert len (meteoPaths) == 1

def paramFile (param):
  """
  create a new .orig file for writing, deleting existing
  :param param: one of "P", "S", etc
  :return:open file
  """
  return open (os.path.join (args.outdir, param + ".orig"), "w")

#create the output files
S_orig = paramFile ("S")
T_orig = paramFile ("T")
H_orig = paramFile ("H")
P_orig = paramFile ("P")
W_orig = paramFile ("W")
D_orig = paramFile ("D")

#empty files, create them and close them
E_orig = paramFile ("E")
R_orig = paramFile ("R")
TETA0_orig = paramFile ("TETA0")
E_orig.close ()
R_orig.close ()
TETA0_orig.close ()

#make S.orig
# #JulianDate	Airmass	na   	na	  na    	na	 na	  seeing
# 2453644.493	1.13	  11730	0.807	0.758	29091	41.4	0.737
inFile = open (seeingPaths [0])
lines = ambient.readJulianDayLines (inFile)
while (lines):
  S_orig.write ("{0}\n".format (len (lines)))
  for line in lines:
    fields = line.split (' ')
    assert len (fields) == 8
    S_orig.write ("{} {}".format (fields [0], fields [7]))
  lines = ambient.readJulianDayLines (inFile)
S_orig.close ()

#make T, H, P and W
# #JulianDate	P         H 	T	         D	           W
# 2453602.307	29.967	  13  57.799999	 15	          40
inFile = open (meteoPaths [0])
lines = ambient.readJulianDayLines (inFile)
while (lines):
  P_orig.write ("{0}\n".format (len (lines)))
  H_orig.write ("{0}\n".format (len (lines)))
  T_orig.write ("{0}\n".format (len (lines)))
  W_orig.write ("{0}\n".format (len (lines)))
  D_orig.write ("{0}\n".format (len (lines)))

  for line in lines:
    fields = line.split (' ')
    assert len (fields) == 6
    date = fields [0]
    P = float (fields [1]) *  24.76520098
    H = fields [2]
    T = (float (fields [3])-32.0)*5.0/9.0
    D = fields [4]
    W = float (fields [5]) * 0.44704

    P_orig.write ("{} {}\n".format (date, P))
    H_orig.write ("{} {}\n".format (date, H))
    T_orig.write ("{} {}\n".format (date, T))
    D_orig.write ("{} {}\n".format (date, D))
    W_orig.write ("{} {}\n".format (date, W))

  lines = ambient.readJulianDayLines (inFile)
P_orig.close ()
H_orig.close ()
T_orig.close ()
W_orig.close ()

