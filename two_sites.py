#!/usr/local/bin/python3

import os
import ambient
import utils

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

parser = ambient.ArgumentParser (description="""
make some plots out of data from 2 sites.
usage: --indirs <sitedir1,startdate1,stopdate1>     <sitedir2,startdate2,stopdate2> --outdir <dir> --min-points <min number of points for pc>
sitedirs are directories containing sets of .orig files
dates are like 1-21-2006
each indir paramter is a comma separated list with no spaces

""")

#parse user's command line arguments
args = parser.parse_args ()

#make sure the output directory exists
utils.makeDir (args.outdir)

#must specify exactly 2 site dirs
assert len (args.indirs) == 2


# 1a select period—1 year, 3mo, etc, user selects
# 1b cycle through each day of period for campanas seeing data. If seeing value 0<S<5, for at least <100> points, record data. Otherwise, label that day ‘no data’
# 1c. record maximum and minimum absolute time of valid data points
# 1d.record valid data points for valid days at lasilla

def getSeeingDays (which):
  """
  :param which: index of indirs
  :return: list of days for seeing values that's got each of its days scanned for a good period of measurements.
  """
  days = ambient.Days (os.path.join (args.indirs [which].path, ambient.parameterBasenames [0] + ".orig"), args.start_date, args.stop_date)
  for day in days.days: day.computeGoodness (0, 5, args.min_points)

  return days

#make lists of all seeing days for each site and subranges of days we are interested in
#see if any of the date ranges is for a leap year so we can pad any non leap year ranges
#with a fake leap day
allSeeing = []
seeing = []
leapYear = False
for i in range (2):
  allSeeing.append (getSeeingDays (i))
  if allSeeing [-1].isLeapYear ():
    leapYear = True

for i in range (2):
  dayList = ambient.dayListFromRange (allSeeing [i], args.indirs [i].startDate, args.indirs [i].stopDate, leapYear=leapYear)
  seeing.append (dayList)


# 1e. calculate PC for each day. Sort from highest to lowest PC
#make a list of (day0, day1, pc)
#day0 is never None but day1 can be None whenever it's lacking a day that site0 has
#pc values > 1.0 mean nonsense, including not enough matching points
correlations = ambient.goodCorrelate (seeing [0], seeing [1], args.min_points)

#sort on pc
correlations = sorted (correlations, key=lambda x: x[2], reverse=not args.worst)
def writeCorr (out, c):
  day0 = c [0]
  day1 = c [1]
  if not day0:
    return

  if not day1 or not day0.good:
    out.write ("{:%m/%d/%y}\n".format (day0.datetime))
  else:
    out.write ("{:%m/%d/%y} {:%m/%d/%y} {:.5f} \
    ({} {}) ({} {}))\n".format (day0.datetime,
                            day1.datetime,
                            c[2],
                            day0.firstGoodMeasurement[0],
                            day1.lastGoodMeasurement[0],
                            day1.firstGoodMeasurement[0],
                            day1.lastGoodMeasurement[0]))


utils.dumpList (os.path.join (args.outdir, "{}-{}.corr").format (args.indirs [0].name, args.indirs [1].name), correlations, writeCorr)

pdfPages = PdfPages (os.path.join (args.outdir, "seeing-{}x{}.pdf").format (args.indirs [0].name, args.indirs [1].name))
matrix = ambient.goodSquareCorrelate (seeing [0], seeing [1], args.min_points)
utils.plot2dArray (matrix, pdfPages,
                  "correlation of \"seeing\" {} &\n {}".format (args.indirs [0], args.indirs [1]),
                  "day", "day", "pearson coefficient")
pdfPages.close ()



# p1. plot days across X axis, and 5 sensors ‘X’ on Y axis, like pc-seeing+all.pdf., except color
#missing data points white to distinguish no data from zero correlation data
#make a 2d color plot of days on the x axis and sensors names on the Y axis and each
#grid point (x,y) is pc seeing vs y at day x

def plotSensorSeeingPC (which):
  days = seeing [which]
  allDays = allSeeing [which]
  siteRange = args.indirs [which]

  yTicks = []
  correlationsArray = []
  for p in ambient.parameterBasenames [1:]:
    otherDays = ambient.Days (os.path.join (siteRange.path, p + ".orig"), args.start_date, args.stop_date)
    otherDays.setGoodness (allDays)
    yTicks.append (p)
    dayList = ambient.dayListFromRange (otherDays, siteRange.startDate, siteRange.stopDate)
    correlations = ambient.goodCorrelate (days, dayList, args.min_points)
    correlationsArray.append (correlations)

  ambient.plotSensorSeeingPC (os.path.join (args.outdir, "pc-seeing+all-{}.pdf").format (siteRange.name),
                              correlationsArray, yTicks, days,
                              args.indirs [which].name)


#plotSensorSeeingPC (0)
plotSensorSeeingPC (1)

# 1f .. select data for P, H, T,W, D (‘X’) at campanas, and at lasilla, for those periods
#
# 1g calculate PC between ‘X’ at Campanas and at lasilla
#
# Make (3X2) plotlets, showing 2 curves in each of 6 panes
#
# 1: S(Campanas), S(lasilla), label day, PC
#
# 2-6 ‘X’ (Campanas),’X’(lasilla)

def plotDay (day, ax, color, linestyle, markerfacecolor, marker):
  measurements = day.goodMeasurements ()
  ax.set_xticklabels([])

  #make x axis absolute
  startSeconds = 0
  ax.scatter ([m[0] - 0 for m in measurements],
           [m[1] for m in measurements],
            edgecolor=color, facecolors="none", marker=marker, s=3)

def setSeeingY (ax):
  ax.set_ylabel ("Seeing")
  ax.set_ylim (*ambient.parameterLimits [0])

def plotSites (basename, correlations):

  assert (args.best_plots <= args.best_count)

  #load in days for every non seeing sensor for both sites in "sites"
  #sites [i][j] is sensor j Days for site i
  #set the goodness of each day to intersect the goodness for the
  #corresponding day in seeing
  sites = []
  for i in range (2):
    sites.append ([])
    for p in ambient.parameterBasenames [1:]:
      sites [-1].append (ambient.Days (os.path.join (args.indirs [i].path, p + ".orig"), args.start_date, args.stop_date))
      sites [-1][-1].setGoodness (allSeeing [i])

  row = 0
  for correlation in correlations:
    if row >= args.best_plots: break
    if correlation [2] >= 1.0: continue

    pdfPages = PdfPages (os.path.join (args.outdir, "{}.{}.pdf".format (basename, row)))
    numPlotted = 0
    fig = None

    def checkPage (last=False):
      nonlocal pdfPages
      nonlocal numPlotted
      nonlocal fig
      numPlotted += 1
      if numPlotted >= 4 or last:
        numPlotted = 0
        fig.subplots_adjust (wspace=0.5, hspace=0.3)
        pdfPages.savefig (fig)
        plt.close (fig)
        if not last:
          fig = plt.figure ()
        else:
          fig = None

    fig = plt.figure ()
    if correlation [1]:
      # Upper Left: Seeing for both sites
      ax = fig.add_subplot (2, 2, numPlotted + 1)
      ax.set_xlabel ("Time")
      setSeeingY (ax)
      ax.set_title ("{:%m/%d/%y},{:%m/%d/%y} pc={:0.3f}".format (correlation [0].datetime, correlation [1].datetime, correlation [2]))
      plotDay (correlation [0], ax, 'red', 'none', 'none', 'o')
      plotDay (correlation [1], ax, 'blue', 'none', 'none', 'v')

    checkPage ()

    assert (len (ambient.parameterBasenames) == len (ambient.parameterLimits))

    for i in range (1, len (ambient.parameterBasenames)):
      p = ambient.parameterBasenames [i]
      limits = ambient.parameterLimits [i]
      otherDays0 = sites [0][i -1]
      otherDays1 = sites [1][i -1]

      otherDay0 = otherDays0.getDay (correlation [0])
      otherDay1 = otherDays1.getDay (correlation [0])

      ax = fig.add_subplot (2, 2,  numPlotted + 1)
      ax.set_ylabel ("{}".format (p))
      ax.set_ylim (*limits)
      ax.set_xlabel ("Time")

      if otherDay0 and otherDay1:
        f0, f1 = [d.goodPeriodFilter () for d in (otherDay0, otherDay1)]
        m0, m1 = otherDay0.matchMeasurements (otherDay1, f0, f1)
        pc = ambient.pearsonFromMeasurements (m0, m1, args.min_points)
        if (pc <= 1):
          ax.set_title ("pc={:0.3f}".format (pc))
        plotDay (otherDay0, ax, 'red', 'none', 'none', 'o')
        plotDay (otherDay1, ax, 'blue', 'none', 'none', 'v')
      checkPage (i == len (ambient.parameterBasenames) - 1)

    pdfPages.close ()

    row += 1


plotSites (("{}-{}.pdf").format (args.indirs [0].name, args.indirs [1].name), correlations)
