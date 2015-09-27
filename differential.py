#!/usr/local/bin/python3

import os
import ambient
import utils
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

parser = ambient.ArgumentParser (description="""
make some plots in a file called differential.pdf showing sensor data in 2 sites vs differences
between seeing values in 2 sites.
usage: --indirs <sitedir1,startdate1,stopdate1>     <sitedir2,startdate2,stopdate2> --outdir <dir> --min-difference=2
sitedirs are directories containing sets of .orig files
dates are like 1-21-2006
each indir paramter is a comma separated list with no spaces
there should be exactly 2 indirs
""")

#parse user's command line arguments
args = parser.parse_args ()

#make sure the output directory exists
utils.makeDir (args.outdir)

#must specify 2 site dirs
assert len (args.indirs) == 2


def collectDayValues (sday0, sday1, vday0, vday1, x0, y0, x1, y1):
  """

  :param sday0: a seeing Day
  :param sday1: matching seeing Day
  :param vday0: a value Day
  :param vday1: matching value Day
  :param x0: list of x values for vday0 (differences)
  :param y0: list of y values for vday0 (values)
  :param x1: list of x values for vday1
  :param y1: list of y values for vday1
  :return:
  """

  #ignore any days where we don't have data for everyone
  if not sday0.good or not sday1.good or not vday0.good or not vday1.good:
    return

  #match up seeing for both days
  any = lambda m: True
  ms0, ms1 = sday0.matchMeasurements (sday1, any, any)
  assert len (ms0) == len (ms1)

  #make list of 3-tuple, (measurement, measurement, abs (difference)
  #for differences above cutoff

  seeingDiffs = []
  for m0, m1 in zip (ms0, ms1):
    diff = m0 [1] - m1 [1]
    if math.fabs (diff) >= args.min_difference:
      seeingDiffs.append ((m0, m1, diff))


  def collectValues (sday, vday, which, x, y):
    """
    :param sday: seeing day
    :param vday: value day
    :param which: 0 or 1, which measurement in the seeingDiff tuple we care about
    :param x: list of x values (differences)
    :param y: list of y values (sensor values)
    :return:
    """

    #match seeing and value measurements for the pair of days
    sm, vm = sday.matchMeasurements (vday, any, any)
    assert len (sm) == len (vm)

    #make a dictionary that matches timestamp values from measurements
    #in sm to corresponding values from measurements in vm.

    timestampToValue = {}
    for m0, m1 in zip (sm, vm):
      timestampToValue [m0 [0]] = m1 [1]

    #for each seeingDiff measurement, if there's a matching value, collect it to x, y
    for t in seeingDiffs:
      value = timestampToValue.get (t[which][0])
      if value != None:
        x.append (t[2])
        y.append (value)

  collectValues (sday0, vday0, 0, x0, y0)
  collectValues (sday1, vday1, 1, x1, y1)


#get 2 lists of days, one list from each site, where the days are matched up
#with "not good" days filling the gaps, where the dates might not have
#matched up
def makeDayLists (paramIndex):
  dayLists = []
  for i in range (2):
    days = ambient.Days (os.path.join (args.indirs [i].path, ambient.parameterBasenames [paramIndex] + ".orig"), args.start_date, args.stop_date)
    dayList = ambient.dayListFromRange (days, args.indirs [i].startDate, args.indirs [i].stopDate)
    dayLists.append (dayList)
  return dayLists

seeingDayLists = makeDayLists (0)
assert len (seeingDayLists [0]) == len (seeingDayLists [1])

def makeScatter (ax, paramIndex):
  valueDayLists = makeDayLists (paramIndex)
  assert len (valueDayLists [0]) == len (valueDayLists [1])
  assert len (seeingDayLists [0]) == len (valueDayLists [0])
  x0 = []
  y0 = []
  x1 = []
  y1 = []
  labels = []
  for i in range (2):
    labels.append (args.indirs [i].name)

  for sday0, sday1, vday0, vday1 in zip (seeingDayLists [0], seeingDayLists [1], valueDayLists [0], valueDayLists [1]):
    collectDayValues (sday0, sday1, vday0, vday1, x0, y0, x1, y1)

  ax.set_title ("{} when Seeing difference is > {:0.3f}\n{:%m/%d/%y}-{:%m/%d/%y}".format (
    ambient.parameterNames [paramIndex],
    args.min_difference,
    args.indirs [0].startDate,
    args.indirs [0].stopDate).title ())

  ax.set_xlabel ("Seeing Difference, {} - {}".format (
    args.indirs [0].name,
    args.indirs [1].name))

  ax.set_ylabel (ambient.parameterNames [paramIndex] + " (" + ambient.parameterUnits [paramIndex] + ")")
  ax.set_ylim (*ambient.parameterLimits [paramIndex])

  p0 = ax.scatter (x0, y0, edgecolor="red", facecolors="none", marker="o", s=3)
  p1 = ax.scatter (x1, y1, edgecolor="blue", facecolors="none", marker="o", s=3)
  ax.legend ((p0, p1), labels)


# T, W, P, H, D for both sites.
pdfPages = PdfPages (os.path.join (args.outdir, "differential.pdf"))

for paramIndex in [1, 4, 3, 2, 8] :
  fig = plt.figure ()
  ax = fig.add_subplot (1, 1, 1)
  makeScatter (ax, paramIndex)
  pdfPages.savefig (fig)
  plt.close (fig)
pdfPages.close ()

