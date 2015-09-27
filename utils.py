import os.path
import re
import datetime
import scipy.stats.stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#string to date converter for argparse
#it returns a datetime if the string was ok, otherwise throws
#an exception that argparse understands
def dateFromString (s):
  try:
    return datetime.datetime.strptime(s, "%m-%d-%Y")
  except ValueError:
    raise argparse.ArgumentTypeError ("Not a valid date: '{0}'.".format (s));


#----- file stuff
def pathFromString (s):
  s = os.path.expanduser (s)
  s = os.path.expandvars (s)
  return s

#create a directory if it doesn't already exist
#this is used for the outputdir option
def makeDir(path):
  try:
    path = pathFromString (path)
    os.makedirs(path)
  except OSError:
    if not os.path.isdir(path):
      raise

#read next line, skipping over comments and blank lines
def readNextLine (file):
  s = file.readline ()
  while s and (s == "\n" or re.match ("#", s)) :
    s = file.readline ()
  return s

def dumpList (path, list, func):
  """
  dump contents of a list to file at path, overwriting
  :param path: path for file to write to
  :param list:
  :param func: func (file, item) to write an item from list to file
  :return:
  """
  file = open (path, "w")
  for item in list:
    func (file, item)
  file.close ()


#write a contourf plot of array2d to path
#z is a 2d array or list indexed like [x][y]
def contourf (z, path, title, xLabel, yLabel, cbarLabel, show=False):
  fig = plt.figure()
  ax = fig.add_subplot(111)
  cs = ax.contourf(z, levels=np.linspace (min (z), max (z), num=255))
  ax.set_title (title.title ())
  ax.set_xlabel (xLabel.title ())
  ax.set_ylabel (yLabel.title ())
  cbar = fig.colorbar(cs)
  cbar.ax.set_ylabel (cbarLabel.title ())
  fig.savefig (path)
  if show:
    plt.show ()
  plt.close (fig)

#write a contourf plot of array2d to path
#z is a 2d array or list indexed like [x][y]
def plot2dArray (z, pdf, title, xLabel, yLabel, cbarLabel, vmin=-1, vmax=1):
  fig = plt.figure()
  ax = fig.add_subplot(111)
  cmap = plt.get_cmap ('jet')
  cmap.set_over ('white')
  cmap.set_under ('white')
  cs = ax.imshow(z, interpolation='nearest',origin='lower', vmin=vmin, vmax=vmax)
  ax.set_title (title.title ())
  ax.set_xlabel (xLabel.title ())
  ax.set_ylabel (yLabel.title ())
  cbar = fig.colorbar(cs)
  cbar.ax.set_ylabel (cbarLabel.title ())
  pdf.savefig (fig)
  plt.close (fig)

def scatterGraph (x, y, path, title, xLabel, yLabel, yticks=None, xticks=None, show=False):
  fig = plt.figure ()
  ax = fig.add_subplot (111)
  ax.set_title (title.title ())
  ax.set_xlabel (xLabel.title ())
  ax.set_ylabel (yLabel.title ())
  if len (x):
    ax.set_xlim (min (x), max (x))
    ax.set_ylim (min (y), max (y))
  if yticks: ax.set_yticks (yticks)
  if xticks: ax.set_yticks (xticks)

  ax.scatter (x, y)
  fig.savefig (path)
  if show:
    plt.show ()
  plt.close (fig)

#class to encapsulate a dir name followed by 2 dates
class SiteRange (object):
  def __init__ (self, argString):
    """
    :param argString: comma separated string from command line, for instance ../lasilla,1-17-2006,10-17-2006
    :return:
    """
    strings = argString.split (',')
    assert len (strings) == 3
    self.path = pathFromString (strings [0])
    self.startDate = dateFromString (strings [1])
    self.stopDate = dateFromString (strings [2])
    assert self.startDate <= self.stopDate
  def __str__(self):
    return "{} {:%m/%d/%y}-{:%m/%d/%y}".format (self.name, self.startDate, self.stopDate)
  @property
  def name (self):
    return os.path.basename (self.path)

