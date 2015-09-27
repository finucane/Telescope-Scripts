#!/usr/local/bin/python3

import os
import ambient
import utils

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

parser = ambient.ArgumentParser (description="""
make some plots out of data from more than 1 site.
this script produces:
sites-seeing.histogram.pdf - a histogram of seeing at multiple sites
sites-temperature.histogram.pdf - histogram of temperature (deg C)
sites-wind-speed.histogram.pdf - histogram of wind (m/s)
sites-wind-direction.histogram.pdf - histogram of direction (degrees)
sites-humidity.histogram.pdf - histogram of humidity (percent)

usage: --indirs <sitedir1,startdate1,stopdate1>     <sitedir2,startdate2,stopdate2> --outdir <dir>
sitedirs are directories containing sets of .orig files
dates are like 1-21-2006
each indir paramter is a comma separated list with no spaces

""")

#parse user's command line arguments
args = parser.parse_args ()

#make sure the output directory exists
utils.makeDir (args.outdir)

#must specify more than 1 site dir
assert len (args.indirs) >= 2

def makeHistogram (paramIndex):

  filename = "sites-" + ambient.parameterNames [paramIndex].translate (str.maketrans (" ", "-")) + ".histogram.pdf"
  filename = filename.lower ()
  title = "Normalized " + ambient.parameterNames [paramIndex]

  #get a list of seeing days, 1 Days object for each site
  sites = []
  for siteRange in args.indirs:
    sites.append (ambient.Days (os.path.join (siteRange.path, ambient.parameterBasenames [paramIndex] + ".orig"), args.start_date, args.stop_date))


  pages = PdfPages (os.path.join (args.outdir, filename))
  fig = plt.figure ()
  ax = fig.add_subplot (1, 1, 1)
  ax.set_color_cycle (['b','r','g','m','c','y'])
  datasets = []
  labels = []

  #for each site, add a list of values to "datasets"
  #also make a corresponding list of labels based on site dir names
  for site, siteRange in zip (sites, [siteRange for siteRange in args.indirs]):
    labels.append (siteRange)
    datasets.append ([])
    dayList = ambient.dayListFromRange (site, siteRange.startDate, siteRange.stopDate)
    for day in dayList:
      if day:
        func = day.periodFilter (args.start_minutes * 60, args.stop_minutes * 60, args.past_midnight)
        filtered = filter (func, day.measurements)
        datasets [-1].extend ([m [1] for m in filtered])

  ax.hist (datasets, bins=args.num_bins, label=labels, normed=1)

  ax.set_title (title.title ())
  ax.set_xlabel (ambient.paramLabel (paramIndex))
  ax.legend ()

  pages.savefig (fig)
  plt.close (fig)
  pages.close ()

# Seeing, Temperature, Wind, Direction, Humidity.
#0, 1, 4, 8, 2
for paramIndex in [0, 1, 4, 8, 2] :
  makeHistogram (paramIndex)
