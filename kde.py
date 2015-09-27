#!/usr/local/bin/python3
#https://github.com/scipy/scipy/blob/master/doc/source/tutorial/stats/plots/kde_plot5.py

import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import ambient
import utils
import numpy as np
import matplotlib.animation as manimation
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages

parser = ambient.ArgumentParser (description="""
kernel density estimation stuff
usage: --indir <sitedir> --outdir <dir> --vars=<var>,<var> --date=<date>
or for a movie:
       --indir <sitedir> --outdir <dir> --vars=<var>,<var> --start-date=<date> --stop-date=<date>

sitedir is a directory containing sets of .orig files
dates are like 1-21-2006
vars is a list of names of any 2 sensors: S, T, H, P, W, E, R, TETA0, D
for example
kde.py --outdir=../scratch --indir=../campanas-all --vars t e --date=5-23-2008
""")

#parse user's command line arguments
args = parser.parse_args ()

#make sure the output directory exists
utils.makeDir (args.outdir)

#get date range
if args.date:
  startDate = stopDate = args.date
else:
  startDate = args.start_date
  stopDate = args.stop_date

#get days objects for each of the 2 variables
#also get min/max values of measurement values
days = []
minmaxes = [None, None]
assert len (args.vars) == 2
for i in range (2):
  days.append (ambient.Days (os.path.join (args.indir, ambient.parameterBasenames [args.vars [i]] + ".orig"), startDate, stopDate))
  for day in days [-1].days:
    minmaxes [i] = day.minMaxValue (minmaxes [i], lambda x:True)

#the rest of the script will just use xmin, xmax, ymin, ymax
xmin = minmaxes [0][0]
xmax = minmaxes [0][1]
ymin = minmaxes [1][0]
ymax = minmaxes [1][1]

#get a pair of lists of day objects for each Days object, matching each day
#gaps will be None's in the 2nd array
matchedDays = ambient.matchDays (days [0].days, days [1].days)


#draw frame, returning True if anything was actually drawn
def draw_frame (fig, im, lines, day0, day1, xmin, xmax, ymin, ymax):
  """
  using set_data and set_array, modify underlying data for an imshow image and plot's lines
  :param fig: figure containing im and lines
  :param im: return value for ax.imshow to shove a matrix into
  :param lines: return value of plot to shove lists of measurements into
  :param m0: measurement array
  :param m1: measurement array matching m1
  :return:
  """

  if not day0 or not day1:
    return False

  #match their measurements
  m0, m1 = day0.matchMeasurements (day1)

  assert not len (m0) or (m0 [0] in day0.measurements)

  #get lists of sensor values from each day
  measurements = []
  measurements.append (np.array ([m [1] for m in m0]))
  measurements.append (np.array ([m [1] for m in m1]))

  X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
  positions = np.vstack([X.ravel(), Y.ravel()])
  values = np.vstack([measurements [0], measurements [1]])
  kernel = stats.gaussian_kde(values)
  Z = np.reshape(kernel.evaluate (positions).T, X.shape)
  im.set_array (np.rot90 (Z))
  lines.set_data (measurements [0], measurements [1])
  return True

#set fig, im, and lines to call draw_frame repeatedly with
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)

im = ax.imshow([[],[]], cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax], aspect='auto')
lines, = ax.plot ([], [], 'k.', markersize=2)
ax.set_title ("Kernel Density Estimation \n {} - {}".format(
  os.path.basename (args.indir),
  startDate.strftime ("%m/%d/%Y")))

ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
ax.set_xlabel (ambient.paramLabel (args.vars [0]))
ax.set_ylabel (ambient.paramLabel (args.vars [1]))

#the script's big switch, either we have 1 date and we make a pdf or its a date range
#and we are making an mp4
if args.date:
  pdf = PdfPages (os.path.join (args.outdir, "kde_{}_{}_{}.pdf".format (
    args.date.strftime ("%m_%d_%Y"),
    ambient.parameterBasenames [args.vars [0]],
    ambient.parameterBasenames [args.vars [1]])))

  draw_frame (fig, im, lines, matchedDays [0][0], matchedDays [1][0], xmin, xmax, ymin, ymax)

  pdf.savefig (fig)
  plt.close (fig)
  pdf.close ()
else:
  assert args.start_date and args.stop_date
  assert len (matchedDays [0])
  FFMpegWriter = manimation.writers['ffmpeg']
  writer = FFMpegWriter (fps=args.fps)
  with writer.saving (fig, os.path.join (args.outdir, "kde_{}_{}_{}_{}.mp4".format (
    startDate.strftime ("%m_%d_%Y"),
    stopDate.strftime ("%m_%d_%Y"),
    ambient.parameterBasenames [args.vars [0]],
    ambient.parameterBasenames [args.vars [1]])), 100):
    for day0, day1, dayCount in zip (matchedDays [0], matchedDays [1], range (len (matchedDays [0]))):
      if draw_frame (fig, im, lines, day0, day1, xmin, xmax, ymin, ymax):
        writer.grab_frame ()
      sys.stdout.write("\r%.1f%%" % (100 * (dayCount / len (matchedDays [0]))))
      sys.stdout.flush()
