#!/usr/local/bin/python3

#fetch_ambient.py

#Y-m-d
#http://archive.eso.org/asm/ambient-server?site=paranal&night=2009-11-01&data=S")
import utils
import ambient

import re
import os.path
import sys
import urllib.request
import http.client
import datetime
import argparse
import dateutil.relativedelta


#make a default string representing some date in the past


#setup an args parser for the script's parameters
parser = ambient.ArgumentParser (description="""Dump ambient sensor values from the eso server to some files.
8 files are created, named S.orig, T.orig, etc, and no processing is done from the original server
data except that the values are sorted chronologically and each day is proceeded by a header line, which is just
the number of samples in the following day, to make it easy to parse the file.
""")

#parse user's command line arguments
args = parser.parse_args ()

#make sure the output directory exists
utils.makeDir (args.outdir)

#fetch contents of a url and append it to the output file 'file'
#on any network error the script fails. urls that correspond
#to no data are ignored
def dump (url, file):
  #add all the lines of data to a list, prepending to reverse the order
  #so its sorted chronologically. if there's not any data, or the data
  #otherwise looks wrong, the function returns before the file
  #is written to.
  lines = []
  response = None
  while True:
    try:
      response = urllib.request.urlopen(url, timeout=180)
      break
    except:
      continue;

  for line in response:
    line = line.decode('utf-8')
    if re.match ("# No data available", line):
      return
    fields = line.split (' ')
    if not len (fields) == 2:
      return
    try:
      f0 = float (fields [0])
      f1 = float (fields [1])
    except:
      return
    lines.insert(0, line);
  file.write ("{0}\n".format (len (lines)))
  file.writelines (lines)


#for each ambient parameter, fetch the date range of data and write to a file
#on any errors throw an exception and abort the script
def fetch (start_date, stop_date, outdir, parameters):
  timedelta = stop_date - start_date
  totalDays = timedelta.days * len (parameters)
  if totalDays == 0: totalDays = len (parameters)
  dayCount = 0
  for p, basename in zip (parameters, ambient.parameterBasenames):
    day = start_date
    f = open (os.path.join (outdir, basename) + ".orig", "w")
    while day <= stop_date:
      url = 'http://archive.eso.org/asm/ambient-server?site={0}&night={1}&data={2}'.format (args.site, day.strftime("%Y-%m-%d"), p)
      sys.stdout.write("\r%.1f%%" % (100 * (dayCount / totalDays)))
      sys.stdout.flush()
      #print (url)
      dump (url, f)
      day += dateutil.relativedelta.relativedelta (days=1)
      dayCount += 1
    f.close ()

assert args.start_date and args.stop_date

fetch (args.start_date, args.stop_date, args.outdir, ambient.parameters)