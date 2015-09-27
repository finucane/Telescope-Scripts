#!/usr/local/bin/python3

import ambient
import utils
import argparse
import os.path
import itertools
import scipy.stats.stats

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from enum import Enum

#----script starts here

parser = ambient.ArgumentParser (description="""generate plots from ESO data""")

#parse user's command line arguments
args = parser.parse_args ()

#make sure the output directory exists
utils.makeDir (args.outdir)

#get list of all seeing days
allSeeingDays = ambient.Days (os.path.join (args.indir, ambient.parameterBasenames [0] + ".orig"), args.start_date, args.stop_date)

#scripts need start/stop dates now
#assert args.start_date

#if start-date/stop-date is specified on command line, only look at a subset of days
if args.start_date:
  seeingDays = ambient.dayListFromRange (allSeeingDays, args.start_date, args.stop_date)
else:
  seeingDays = allSeeingDays.days


def seeingCorrelations (seeingDays, otherDays):
  """
  :param seeingDays: a list of Day objects for seeing
  :param otherDays: a list of Day objects for some sensor
  :return: correlation list of 3-tuples of pc's
  """
  correlations = ambient.crossCorrelate (seeingDays, otherDays, args.start_minutes*60, args.stop_minutes*60, args.past_midnight, args.min_points)

  #this is a 3-tuple (day, otherDay, pc), day is in seeingDays, otherDay is in otherDays

  #remove rows with nonsense pc values
  correlations = [t for t in correlations if t[2] <= 1.0]

  #filter just pcs in the range we want
  correlations = [t for t in correlations if t [2] >= args.start_pc and t [2] <= args.stop_pc]

  return correlations

#make a 2d color plot of days on the x axis and sensors names on the Y axis and each
#grid point (x,y) is pc seeing vs y at day x
yTicks = []
correlationsArray = []
for p in ambient.parameterBasenames [1:]:
  otherDays = ambient.Days (os.path.join (args.indir, p + ".orig"), args.start_date, args.stop_date, allSeeingDays)
  yTicks.append (p)
  correlations = ambient.crossCorrelate (seeingDays, otherDays.days, args.start_minutes*60, args.stop_minutes*60, args.past_midnight, args.min_points)
  correlationsArray.append (correlations)

ambient.plotSensorSeeingPC (os.path.join (args.outdir, "pc-seeing+all.pdf"), correlationsArray, yTicks, seeingDays, os.path.basename (args.indir))

#make histogram of seeing for all seconds of days and
#seeing for start/stop seconds of all days
pdfPages = PdfPages (os.path.join (args.outdir, "seeing.histogram.pdf"))
fig = plt.figure ()
ax = fig.add_subplot (1, 1, 1)

allSeeing = []
periodSeeing = []
for day in seeingDays:
  values = [m [1] for m in day.measurements]
  allSeeing.extend (values)

  func = day.periodFilter (args.start_minutes * 60, args.stop_minutes * 60, args.past_midnight)
  filtered = filter (func, day.measurements)
  periodSeeing.extend ([m [1] for m in filtered])

ax.hist ([allSeeing, periodSeeing], bins=args.num_bins, label=['All Day','Time Range'])

ax.set_title ("All seeing values and seeing values for daily time range".title ())
ax.set_xlabel ("Seeing")
ax.set_ylabel ("Count".title ())
ax.legend ()

pdfPages.savefig (fig)
plt.close (fig)
pdfPages.close ()

#get list of (day, day, pc)
correlations = ambient.selfCorrelate (seeingDays, args.start_minutes*60, args.stop_minutes*60, args.past_midnight, args.min_points)
#sort on pc
correlations = sorted (correlations, key=lambda x: x[2], reverse=not args.worst)

utils.dumpList (os.path.join (args.outdir, "SxS.corr"), correlations,
               lambda out, item: out.write ("{:%m/%d/%y} {:%m/%d/%y} {:.5f} ({:d}, {:d}, {:d})\n".format (item[0].datetime, item[1].datetime, item [2],
                                                                                        item[3][0],
                                                                                        item[3][1],
                                                                                        item[3][2])))

#remove rows with nonsense pc values
#leave them in and they'll be colored white
#correlations = [t for t in correlations if t[2] <= 1.0]

#make a matrix instead of a list so we can do a matrix plot of everything
corrMatrix = ambient.makeCorrelationMatrix (seeingDays, correlations)
pdfPages = PdfPages (os.path.join (args.outdir, "SxS.pdf"))

utils.plot2dArray (corrMatrix, pdfPages,
                  "self correlation of \"seeing\" at {}\n{:%m/%d/%Y} - {:%m/%d/%Y}".
                  format (os.path.basename (args.indir), seeingDays[0].datetime, seeingDays[-1].datetime),
                  "day", "day", "pearson coefficient")
pdfPages.close ()


#filter just pcs in the range we want
correlations = [t for t in correlations if t [2] >= args.start_pc and t [2] <= args.stop_pc]

#take subset of best (or worst) pcs
if args.best_count > 0 and args.best_count < len (correlations):
  correlations = correlations [0:args.best_count + 1]

#add a rank column
correlations = [(index, t[0], t[1], t[2]) for index, t in enumerate (correlations)]


utils.dumpList (os.path.join (args.outdir, "SxS.subset.corr"), correlations,
               lambda out, item: out.write ("{} {:%m/%d/%y} {:%m/%d/%y} {:.5f}\n".format (item [0], item[1].datetime, item[2].datetime, item [3])))

#slice off lists of days from tuples in the correlations list
days0 = [t[1] for t in correlations]
days1 = [t[2] for t in correlations]

#make subplots of SxS


#first get min,max values for seeing values and timestamps (the Y and X axises)
#only look at measurements in the start/stop_minutes range

numPlots = min (args.num_subplots, len (correlations))
numPlotsPerPage = min (numPlots, args.num_rows * args.num_cols)

minMaxValue = None
minMaxUtfSeconds = (args.start_minutes*60, args.stop_minutes*60)

for i in range (numPlots):
  day0 = correlations [i][1]
  day1 = correlations [i][2]
  func = day0.periodFilter (args.start_minutes * 60, args.stop_minutes * 60, args.past_midnight)
  minMaxValue = day0.minMaxValue (minMaxValue, func)

  func = day1.periodFilter (args.start_minutes * 60, args.stop_minutes * 60, args.past_midnight)
  minMaxValue = day1.minMaxValue (minMaxValue, func)

#chop off measurements out of range and then correct the x-axis to be zeroed at the period start
def plotDay (day, ax, color, linestyle, markerfacecolor, marker):
  func = day.periodFilter (args.start_minutes * 60, args.stop_minutes * 60, args.past_midnight)
  measurements = list (filter (func, day.measurements))
  if args.past_midnight:
    startSeconds = 0
  else:
    startSeconds = measurements [0][0]
  ax.set_xticklabels([])

  #make x axis absolute
  startSeconds = 0
  ax.scatter ([m[0] - startSeconds for m in measurements],
           [m[1] for m in measurements],
            edgecolor=color, facecolors="none", marker=marker, s=3)

def setSeeingY (ax):
  ax.set_ylabel ("Seeing")
  ax.set_ylim (*ambient.parameterLimits [0])

pdfPages = PdfPages (os.path.join (args.outdir, "SxS.best.scatter.pdf"))

#plot each subplot
numPlotsWritten = 0
for i in range (numPlots):
  day0 = correlations [i][1]
  day1 = correlations [i][2]

  #also make a new figure when we get to the start of a new page
  if (i % numPlotsPerPage) == 0:
    fig = plt.figure ()
  ax = fig.add_subplot (args.num_rows, args.num_cols, (i % numPlotsPerPage) + 1)

  #label the first subplot more than the rest
  if (i % numPlotsPerPage) == 0:
    ax.set_xlabel ("Seconds")
    setSeeingY (ax)
  else:
    ax.set_yticklabels([])

  #set subplot's title and x,y axis ranges
  ax.set_title ("{:%m/%d/%y},{:%m/%d/%y}".format (day0.datetime, day1.datetime))
  ax.set_ylim (*ambient.parameterLimits [0])
  ax.set_xlim (minMaxUtfSeconds [0], minMaxUtfSeconds [1])

  plotDay (day0, ax, 'red', 'none', 'none', 'o')
  plotDay (day1, ax, 'blue', 'none', 'none', 'v')

  #close off the page if we have just written some multiple of numPlotsPerPage, or if
  #we just wrote the last subplot
  numPlotsWritten += 1
  if (numPlotsWritten == numPlotsPerPage or i == numPlots - 1):
    numPlotsWritten = 0
    fig.subplots_adjust(hspace=0.5)
    pdfPages.savefig (fig)
    plt.close (fig)
pdfPages.close ()

#make Big Matrix
def makeBigMatrix (correlations):
  m = []
  for t in correlations:
    row = []
    row.append (t [0]) #rank
    row.append (t [1]) #day of p1
    row.append (t [2]) #day of p2
    row.append (t [3]) #self pc
    m.append (row)
  return m

class Axis (Enum):
  x = 0
  y = 1
  both = 2

#make a scatter graph w/ pc values on the y axis
def scatterPC (x, y, path, title, xLabel, yLabel, pcAxis):
  fig = plt.figure ()
  ax = fig.add_subplot (111)
  ax.set_title (title.title ())
  ax.set_xlabel (xLabel.title ())
  ax.set_ylabel (yLabel.title ())

  if pcAxis == Axis.both or pcAxis == Axis.x:
    ax.set_xlim (-1, 1)
    ax.set_xticks (np.arange (-1, 1, 0.25))
  elif len (x):
    ax.set_xlim (min (x), max (x))

  if pcAxis == Axis.both or pcAxis == Axis.y:
    ax.set_ylim (-1, 1)
    ax.set_yticks (np.arange (-1, 1, 0.25))
  elif len (y):
    ax.set_ylim (min (y), max (y))


  ax.scatter (x, y)
  fig.savefig (path)
  plt.close (fig)

#do seeing x all correlations
for index, days in enumerate ([days0, days1]):
  bigMatrix = makeBigMatrix (correlations)
  for p in ambient.parameterBasenames [1:]:
    otherDays = ambient.Days (os.path.join (args.indir, p + ".orig"), args.start_date, args.stop_date, allSeeingDays)

    #for each day in the seeing correlations tuple, compute the pc seeing vs a non-seeing sensor
    for row, day in enumerate (days):
      otherDay = otherDays.getDay (day)
      pc = (0,(0,0,0))
      if otherDay:
        f0, f1 = [d.periodFilter (args.start_minutes * 60, args.stop_minutes * 60, args.past_midnight) for d in (day, otherDay)]
        m0, m1 = day.matchMeasurements (otherDay, f0, f1)
        pc = (ambient.pearsonFromMeasurements (m0, m1, args.min_points), (len (day.measurements), len (otherDay.measurements), len (m0)))
      bigMatrix [row].append (pc)

    #strip out nonsense correlations
    strippedMatrix = [t for t in bigMatrix if t[-1][0] <= 1.0]

    scatterPC ([t [0] for t in strippedMatrix],
                        [t [-1][0] for t in strippedMatrix],
                        os.path.join (args.outdir, "{}x{}.{}.png".format (p, ambient.parameterBasenames [0], index)),
                        "{} vs {} ({})".format (ambient.parameterBasenames [0], p, index),
                        "rank",
                        "pearson's correlation coefficient",
                        Axis.y)

  utils.dumpList(os.path.join (args.outdir, "{}xAll.{}.corr".format(ambient.parameterBasenames[0], index)), bigMatrix,
                 lambda out, r: out.write("{} {:%m/%d/%y} {:%m/%d/%y} {:.5f} {:.5f} ({:d},{:d},{:d}) {:.5f} ({:d},{:d},{:d}) {:.5f} ({:d},{:d},{:d}) {:.5f} ({:d},{:d},{:d}) {:.5f} ({:d},{:d},{:d}) {:.5f} ({:d},{:d},{:d}) {:.5f} ({:d},{:d},{:d})\n".format
                                          (r[0],
                                           r[1].datetime,
                                           r[2].datetime,
                                           r[3],
                                           r[4][0],r[4][1][0],r[4][1][1],r[4][1][2],
                                           r[5][0],r[5][1][0],r[5][1][1],r[5][1][2],
                                           r[6][0],r[6][1][0],r[6][1][1],r[6][1][2],
                                           r[7][0],r[7][1][0],r[7][1][1],r[7][1][2],
                                           r[8][0],r[8][1][0],r[8][1][1],r[8][1][2],
                                           r[9][0],r[9][1][0],r[9][1][1],r[9][1][2],
                                           r[10][0],r[10][1][0],r[10][1][1],r[10][1][2])))



#do other x other correlations
bigMatrix = makeBigMatrix (correlations)
for p in ambient.parameterBasenames [1:]:
  otherDays = ambient.Days (os.path.join (args.indir, p + ".orig"), args.start_date, args.stop_date, allSeeingDays)

  #for each pair of days in the seeing correlations tuple, compute the pc of a non-seeing vs itself
  row = 0
  for day0, day1 in zip (days0, days1):
    #find matching pair of days in the otherDays list
    otherDay0 = otherDays.getDay (day0)
    otherDay1 = otherDays.getDay (day1)

    pc = (0,(0,0,0))
    if otherDay0 and otherDay1:
      f0, f1 = [d.periodFilter (args.start_minutes * 60, args.stop_minutes * 60, args.past_midnight) for d in (otherDay0, otherDay1)]
      m0, m1 = otherDay0.matchMeasurements (otherDay1, f0, f1)
      pc = (ambient.pearsonFromMeasurements (m0, m1, args.min_points), (len (otherDay0.measurements), len (otherDay1.measurements), len (m0)))
    bigMatrix [row].append (pc)
    row += 1

  #strip out nonsense correlations
  strippedMatrix = [t for t in bigMatrix if t[-1][0] <= 1.0]

  scatterPC ([t [3] for t in strippedMatrix],
                      [t [-1][0] for t in strippedMatrix],
                      os.path.join (args.outdir, "{}{}x{}{}.png".format (p, p, ambient.parameterBasenames[0], ambient.parameterBasenames [0])),
                      "{} vs {}".format (p, ambient.parameterBasenames [0]),
                      "{} self correlation".format (ambient.parameterBasenames [0]),
                      "{} self correlation".format (p),
                      Axis.both)

utils.dumpList(os.path.join (args.outdir, "{}{}xAllAll.corr".format (ambient.parameterBasenames[0], ambient.parameterBasenames[0])), bigMatrix,
                lambda out, r: out.write("{} {:%m/%d/%y} {:%m/%d/%y} {:.5f} {:.5f} ({:d},{:d},{:d}) {:.5f} ({:d},{:d},{:d}) {:.5f} ({:d},{:d},{:d}) {:.5f} ({:d},{:d},{:d}) {:.5f} ({:d},{:d},{:d}) {:.5f} ({:d},{:d},{:d})\n".format
                                          (r[0],
                                           r[1].datetime,
                                           r[2].datetime,
                                           r[3],
                                           r[4][0],r[4][1][0],r[4][1][1],r[4][1][2],
                                           r[5][0],r[5][1][0],r[5][1][1],r[5][1][2],
                                           r[6][0],r[6][1][0],r[6][1][1],r[6][1][2],
                                           r[7][0],r[7][1][0],r[7][1][1],r[7][1][2],
                                           r[8][0],r[8][1][0],r[8][1][1],r[8][1][2],
                                           r[9][0],r[9][1][0],r[9][1][1],r[9][1][2],
                                           r[10][0],r[10][1][0],r[10][1][1],r[10][1][2])))



#For each of the top <10> pairs of seeing periods, display 7 charts, 1 for each environmental variable, each comprising 4 subplots,
# each subplot having two variables represented by unfilled circles and unfilled triangles.
# The four subplots are
#
# Upper Left: Seeing for periods 1 and 2.
# UR: Environmental variable for periods 1 and 2.
# LL: Dual ordinate, Seeing and Environmental variable, period 1.
# LR: Dual ordinate, Seeing and Environmental variable, period 2.
# x-axis should be p1 or p2; left y-axis alike, right y-axis alike on each chart, scaled appropriately to each variable.

def quadPlot (days0, days1, basename, correlations):

  assert (args.best_plots <= args.best_count)

  row = 0
  for day0, day1 in zip (days0, days1):
    if row >= args.best_plots: break

    pdfPages = PdfPages (os.path.join (args.outdir, "{}.{}.pdf".format (basename, row)))

    assert (len (ambient.parameterBasenames) == len (ambient.parameterLimits))

    #make 7 pages in the pdf file, one page per sensor
    for paramIndex in range (len (ambient.parameterBasenames) - 1):
      p = ambient.parameterBasenames [paramIndex + 1]
      limits = ambient.parameterLimits [paramIndex +1]
      label = ambient.parameterNames [paramIndex + 1]
      units = ambient.parameterUnits [paramIndex + 1]

    #for p, limits, label in zip (ambient.parameterBasenames [1:], ambient.parameterLimits [1:], amb):
      otherDays = ambient.Days (os.path.join (args.indir, p + ".orig"), args.start_date, args.stop_date, allSeeingDays)

      #find matching pair of days in the otherDays list
      otherDay0 = otherDays.getDay (day0)
      otherDay1 = otherDays.getDay (day1)

      fig = plt.figure ()

      #some .orig files are empty
      if otherDay0 and otherDay1:

        # Upper Left: Seeing for periods 1 and 2.
        ax = fig.add_subplot (2, 2, 1)
        ax.set_xlabel ("Seconds")
        setSeeingY (ax)
        ax.set_title ("{:%m/%d/%y},{:%m/%d/%y} pc={:0.3f}".format (day0.datetime, day1.datetime, correlations [row][3]))
        plotDay (day0, ax, 'red', 'none', 'none', 'o')
        plotDay (day1, ax, 'blue', 'none', 'none', 'v')

        ylabel = label + " (" + units + ")"
        #UR: Environmental variable for periods 1 and 2.
        ax = fig.add_subplot (2, 2, 2)
        ax.set_ylabel (ylabel)
        ax.set_ylim (*limits)
        ax.set_xlabel ("Seconds")
        plotDay (otherDay0, ax, 'red', 'none', 'none', 'o')
        plotDay (otherDay1, ax, 'blue', 'none', 'none', 'v')

        #LL: Dual ordinate, Seeing and Environmental variable, period 1.
        ax = fig.add_subplot (2, 2, 3)
        ax2 = ax.twinx ()
        ax.set_xlabel ("Seconds")
        setSeeingY (ax)
        ax.yaxis.label.set_color('red')
        plotDay (day0, ax, 'red', 'none', 'none', 'o')

        ax2.set_ylim (*limits)
        ax2.set_ylabel (ylabel)
        ax2.yaxis.label.set_color('blue')
        plotDay (otherDay0, ax2, 'blue', 'none', 'none', 'v')

        #LR: Dual ordinate, Seeing and Environmental variable, period 2.
        ax = fig.add_subplot (2, 2, 4)
        ax2 = ax.twinx ()
        ax.set_xlabel ("Seconds")
        setSeeingY (ax)
        ax.yaxis.label.set_color('red')
        plotDay (day1, ax, 'red', 'none', 'none', 'o')

        ax2.set_ylim (*limits)
        ax2.set_ylabel (ylabel)
        ax2.yaxis.label.set_color('blue')
        plotDay (otherDay1, ax2, 'blue', 'none', 'none', 'v')

        fig.subplots_adjust (wspace=0.5, hspace=0.3)
      #close off the page
      pdfPages.savefig (fig)
      plt.close (fig)

    row += 1
    pdfPages.close ()


quadPlot (days0, days1, "seeing+sensors", correlations)


# Step 4. Concentrate on four parameters, Wind Speed, wind direction, humidity, and temperature, “X”. For each “X"
# 4a. For each day Calculate and Rank order PC between Seeing on that day and “X” on that day.
# 4b. Make 4 plots of the top 10 (or so) highest ranking days
# 	1. Seeing vs time, “X” vs time, for that day
# 	2-4. Seeing vs time, not-“X” vs time, for that day
# 4c. Repeat for each “X"

#parameters = ['S', 'T', 'H', 'P', 'W', 'E', 'R', 'TETA0', 'w']

#make list of indexes of the sensors we're interested in
interestingIndexes = [ambient.parameters.index (x) for x in ['W', 'w', 'H', 'T']]

#make list of "Days" objects for each sensor
otherDaysList = [ambient.Days (os.path.join (args.indir, ambient.parameterBasenames [i] + ".orig"), args.start_date, args.stop_date, allSeeingDays) for i in interestingIndexes]

for i in range (len (interestingIndexes)):
  correlateIndex = interestingIndexes [i]

  #make pc correlations between seeing and the next interesting sensor
  otherDays = otherDaysList [i]
  correlations = seeingCorrelations (seeingDays, otherDays.days)

  #sort by pc
  correlations = sorted (correlations, key=lambda x: x[2], reverse=not args.worst)

  pdfPages = PdfPages (os.path.join (args.outdir, "seeing+{}.pdf".format (ambient.parameterBasenames [correlateIndex])))

  for row in range (min (args.best_plots, len (correlations))):

    assert (len (ambient.parameterBasenames) == len (ambient.parameterLimits))

    fig = plt.figure ()
    seeingDay = correlations [row] [0]
    pc = correlations [row] [2]
    fig.suptitle ("{:%m/%d/%y} Seeing and {}, pc={:0.3f}".format (seeingDay.datetime, ambient.parameterBasenames [correlateIndex], pc))

    for j in range (len (interestingIndexes)):
      sensorIndex = interestingIndexes [j]

      otherDay = otherDaysList[j].getDay (seeingDay)

      #some .orig files are empty, also the data sets don't all have the same days
      if not otherDay: continue

      #Seeing and Environmental variable

      ax = fig.add_subplot (2, 2, j + 1)
      ax2 = ax.twinx ()
      ax.set_xlabel ("Time")
      setSeeingY (ax)
      ax.yaxis.label.set_color('red')
      plotDay (seeingDay, ax, 'red', 'none', 'none', 'o')

      ax2.set_ylim (*ambient.parameterLimits [sensorIndex])
      ax2.set_ylabel ("{}".format (ambient.parameterBasenames [sensorIndex]))
      ax2.yaxis.label.set_color('blue')

      plotDay (otherDay, ax2, 'blue', 'none', 'none', 'v')

    fig.subplots_adjust (wspace=0.5, hspace=0.3)
    #close off the page
    pdfPages.savefig (fig)
    plt.close (fig)

  row += 1
  pdfPages.close ()

