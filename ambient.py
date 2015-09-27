import utils

import os
import math
import argparse
import datetime
import itertools
import scipy.stats.stats
import astropy.time
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

#list of ambient parameters, these are used in the server's url syntax and also
#for output filenames. 'S' has to be first since seeing is special case
parameters = ['S', 'T', 'H', 'P', 'W', 'E', 'R', 'TETA0', 'w']
parameterBasenames = ['S', 'T', 'H', 'P', 'W', 'E', 'R', 'TETA0', 'D']
parameterLimits = [(0, 3), (-5, 20), (0, 100), (735, 755), (0, 25), (0, 0.04), (0, 4), (0, 3), (0, 360)]
parameterNames = ["Seeing", "Temperature", "Humidity", "Atmospheric Pressure", "Wind Speed", "E", "R", "TETA0", "Wind Direction"]
parameterUnits = ["", "C", "%", "hPa", "m/s", "", "", "", "Degrees"]

#function to make axis label string
def paramLabel (paramIndex):
  if len (parameterUnits [paramIndex]) > 0:
    return "{} ({})".format (parameterNames [paramIndex], parameterUnits [paramIndex])
  else:
    return "{}".format (parameterNames [paramIndex])


#common arguments across tools
class ArgumentParser (argparse.ArgumentParser):
  def __init__(self, description):
    argparse.ArgumentParser.__init__ (self, description)
    self.sites = ['paranal', 'lasilla']
    self.zeros = ['midnight-ut', 'first-seeing', 'sunrise']

    site = self.sites [0]

    self.add_argument("--start-date", default=None, type=utils.dateFromString,
                    help="first date of data for instance 3-26-2015")
    self.add_argument("--stop-date", default=None, type=utils.dateFromString,
                    help="last date of data for instance 3-26-2015")
    self.add_argument("--date", default=None, type=utils.dateFromString,
                    help="first date of data for instance 3-26-2015")
    self.add_argument ("--site", default=site, choices=self.sites,
                     help="which telescope site, default is {}".format (site))
    self.add_argument ("--otherindir", default=None,
                     help="another input dir for correlations between sites, default is none")
    self.add_argument ("--outdir", default=os.getcwd(), type=utils.pathFromString,
                     help="directory to write files to, default is current directory.")
    self.add_argument ("--indirs", nargs='+', default=None, type=utils.SiteRange,
                     help="list of dir,start-date,stop-date tuples")
    self.add_argument ("--indir", default=os.getcwd(), type=utils.pathFromString,
                     help="directory to read files from, default is current directory.")
    self.add_argument ("--start-minutes", default=0, type=float,
                     help="start of timestamp range, default is 0")
    self.add_argument ("--stop-minutes", default=1440, type=float,
                     help="stop of timestamp range, default is 1440")
    self.add_argument ("--best-count", default=100, type=int,
                     help="how many best items out of a set to pick, default is 100, 0 means all")
    self.add_argument ("--past-midnight", action="store_true",
                     help="turn past midnight on, default is to measure periods from first timestamps in a set")
    self.add_argument ("--worst", action="store_true",
                     help="turn on reversal the sense of best-count")
    self.add_argument ("--start-pc", default=-1, type=float,
                     help="start of pc range, default is -1")
    self.add_argument ("--stop-pc", default=1, type=float,
                     help="stop of pc range, default is 1")
    self.add_argument ("--min-points", default=100, type=int,
                     help="minimum number of points acceptable for pc computations, default is 100")
    self.add_argument ("--min-difference", default=0.0, type=float,
                     help="minimum absolute difference value for differential stuff")
    self.add_argument ("--num-rows", default=3, type=int,
                     help="number of rows per page of subplots, default is 3")
    self.add_argument ("--num-cols", default=2, type=int,
                     help="number of columns per page of subplots, default is 2")
    self.add_argument ("--num-subplots", default=12, type=int,
                     help="number of subplots per pdf file, default is 12")
    self.add_argument ("--best-plots", default=10, type=int,
                     help="number of best data sets to plot default is 10, should be <= --best-count")
    self.add_argument ("--num-bins", default=20, type=int,
                     help="number of bins for histogram, default is 20")
    self.add_argument ("--vars", nargs='+', default=None, type=lambda s: parameterBasenames.index (s.upper ()))
    self.add_argument ("--fps", default=10, type=float,
                       help="frames per second, default is 10")
    self.add_argument ("--bitrate", default=64000, type=int,
                       help="bitrate, default is 64000")


#immutable data type that represents a date without knowing
#its year. it's good for comparisons and hashing and as well
#as knowing a month number and a day number (1-based)

class DayOfYear (tuple):
  def __new__ (cls, month, day):
    return super (DayOfYear, cls).__new__(cls, (month, day))
  def __hash__ (self):
    assert self.month > 0
    assert self.day > 0
    return self.month * 31 + self.day
  def __lt__ (self, other):
    return self.__hash__ () < other.__hash__ ()
  def __le__ (self, other):
    return self.__hash__ () <= other.__hash__ ()
  def __gt__ (self, other):
    return self.__hash__ () > other.__hash__ ()
  def __ge__ (self, other):
    return self.__hash__ () >= other.__hash__ ()
  def __eq__ (self, other):
    return self.__hash__ () == other.__hash__ ()
  def __ne__ (self, other):
    return self.__hash__ () != other.__hash__ ()
  def __str__(self):
    return "{}-{}".format (self.month, self.day)
  @property
  def month (self):
    return self [0]
  @property
  def day (self):
    return self [1]


class Day (object):
  leapDay = DayOfYear (2, 29)

  def __init__ (self, date=None, dayOfYear=None):
    self.measurements = []
    self.cursor = 0
    self.periodZeroUtfSeconds = 0
    self.otherDayWasFound = False #for testing
    self.good = False
    self.goodRange = [-1, -1]
    self.year = None

    self.julianDay = None
    #to match days ignoring year we care about month/day
    #these are 1-based
    self.dayOfYear = None
    self.datetime = None

    #if there was a date specified it's going to be an empty day
    #and we have to set the date
    if date:
      dt = datetime.datetime.combine (date, datetime.time.min)
      assert dt.date () == date
      t = astropy.time.Time (dt, format="datetime", scale='utc')
      julianDay = math.modf (t.jd1) [1]
      self.computeDate (julianDay)

      #mother of all kludges, this happens because of the half-day julian day thing
      #since i hate scripts this is good enough
      if self.datetime.date () != date:
        julianDay += 1
      self.computeDate (julianDay)
      assert self.datetime.date () == date

      #we put this back in, the code is all broken to hell because it's supposed
      #to support matching dates but also matching day-of-years w/out reference
      #to actual years, depending on what scripts are being used
      self.julianDay = julianDay
    elif dayOfYear:
      #this is an empty leap day. datetime cannot be set
      #leave this all broken for now to see if anyone uses datetime
      self.dayOfYear = dayOfYear


  #do all the date computation stuff from julianDay so we can rely on it
  #and not have to keep doing it. this thing should be called as soon
  #as we set julianDay
  def computeDate (self, julianDay):
    self.datetime = astropy.time.Time (julianDay, format="jd").datetime
    date = self.datetime.date ()
    self.dayOfYear = DayOfYear (date.month, date.day)
    self.year = date.year
  #true if day is in a leap year
  def isLeapYear (self):
    assert self.datetime
    return self.year % 4 == 0 and self.year % 100 != 0

  def addMeasurement (self, julianDay, value):
    assert julianDay > 0

    #a Day with a measurement is good until we know better
    self.good = True
    fraction = math.modf (julianDay)
    utfSeconds = (fraction [0] - .5) * 24.0 * 60.0 * 60.0
    #can be -1 if the julian day really was the day before
    self.measurements.append ((utfSeconds, value))

    #set periodZero to the first timestamp
    #also set the date for this day now that we have a timestamp to get the day from
    if len (self.measurements) == 1:
      periodZeroUtfSeconds = self.measurements [0][0]
      self.goodRange [0] = 0
      julianDay = fraction [1]

      #put self.julianDay back in
      self.julianDay = julianDay
      self.computeDate (julianDay)
    self.goodRange [1] = len (self.measurements) - 1

  #go from both ends of the measurements list past
  #points out of bounds on either side, recomputing goodRange
  #and setting "good" if the final goodRange contains at least
  #"minPoints" points
  #goodRange can contain out of bounds points. the idea here is
  #that in the start and end of a day things aren't good, but
  #once the day gets going even with some badness from time to time
  #it's all good
  def setGoodRange (self, func):
    self.goodRange [0] = 0

    self.goodRange [0] = 0
    while self.goodRange [0] < len (self.measurements):
      if func (self.measurements [self.goodRange [0]]):
        break
      self.goodRange [0] += 1

    self.goodRange [1] = len (self.measurements) - 1
    while self.goodRange [1] > self.goodRange [0]:
      if func (self.measurements [self.goodRange [1]]):
        break
      self.goodRange [1] -= 1

  def computeGoodness (self, minValue, maxValue, minPoints):
    assert minPoints > 0
    assert minValue <= maxValue

    def func (m):
      return m [1] >= minValue and m [1] <= maxValue

    self.setGoodRange (func)
    assert self.goodRange [0] <= self.goodRange [1]

    self.good = False
    numGood = 0
    i = self.goodRange [0]
    while i >= 0 and i < len (self.measurements) and i < self.goodRange [1]:
      if func (self.measurements [i]):
        numGood += 1
      i += 1

    self.good = numGood >= minPoints
    assert not self.good or (self.goodRange [0] >= 0 and self.goodRange [0] < len (self.measurements))
    assert not self.good or (self.goodRange [1] >= 0 and self.goodRange [1] < len (self.measurements))

  #set the goodness of self to intersect the goodness of day
  #this means that self's good measurements are the ones that
  #lie within the same time range as that of "day"
  #there are no minPoint restrictions
  def setGoodness (self, day):

    if not day.good:
      self.good = False
      return

    firstTime = day.firstGoodMeasurement [0]
    lastTime = day.lastGoodMeasurement [0]

    def func (m):
       return m [0] >= firstTime and m [0] <= lastTime

    self.setGoodRange (func)
    self.good = self.goodRange [1] >= self.goodRange [0]


  #return list of good measurements, or the empty list if the day is bad
  def goodMeasurements (self):
    if not self.good:
      return []
    return self.measurements [self.goodRange[0] : self.goodRange [1] + 1]


  @property
  def firstGoodMeasurement (self):
    return self.measurements [self.goodRange [0]]
  @property
  def lastGoodMeasurement (self):
    return self.measurements [self.goodRange [1]]

  #return a filter that will accept measurements within the good range, if any
  def goodPeriodFilter (self):
    if not self.good:
      return lambda m: False
    startTime = self.measurements [self.goodRange [0]][0]
    stopTime = self.measurements [self.goodRange [1]][0]
    def func (m):
      return m[0] >= startTime and m[0] <= stopTime
    return func

  def periodFilter (self, startSeconds, stopSeconds, pastMidnight):
    if not self.good:
      return lambda m:False

    start = 0
    if pastMidnight:
      start = 0
    else:
      start = self.periodZeroUtfSeconds
      if not self.otherDayWasFound:
        assert self.otherDayWasFound
    def func (m):
      return m[0] >= start + startSeconds and m[0] <= start + stopSeconds
    return func

  #return 2 lists of measurements that match, based on timestamp, the measurements from self
  #with measurements from an "other" day. the 1st list returned has measurements from self.

  def matchMeasurements (self, otherDay, func=lambda m: True, otherFunc= lambda m: True):

    #throw out bad points
    filtered0 = [m for m in self.measurements if func (m)]
    filtered1 = [m for m in otherDay.measurements if otherFunc (m)]

    #determine which list to match from and which to match to.
    #the from list is the list with fewer points, or the list with the earliest
    #measurement in case of tie. if that also ties then it doesn't matter
    #which list comes first.

    #true if list of measurements ms0's first measurement is before that of ms1
    def earlier (ms0, ms1):
      assert len (ms0) == len (ms1)
      for m0, m1 in zip (ms0, ms1):
        if m0 [0] < m1 [0]:
          return True
        if m0 [0] > m1 [0]:
          return False
      return True #lists happened to have equal times at each measurement

    if len (filtered0) < len (filtered1):
      fromList = filtered0
      toList = filtered1
    elif len (filtered0) > len (filtered1) or not earlier (filtered0, filtered1):
      fromList = filtered1
      toList = filtered0
    else:
      fromList = filtered0
      toList = filtered1

    #return lists so that 1st one's measurements is from self
    def answer (m0, m1):
      if fromList == filtered0:
        return m0, m1
      else:
        return m1, m0

    #get rid of the corner case of fewer than 3 points
    if len (fromList) < 3:
      return answer (fromList, toList [:len (fromList)])


    #determine the boundaries between the points on fromList. these are just times midway between
    #each consecutive measurement

    #start with a placeholder. we'll mirror the first boundary to be backwards from the second
    #across the first point
    boundaries = [0]
    for i in range (len (fromList) - 1):
      m0 = fromList [i]
      m1 = fromList [i + 1]
      assert m0 [0] <= m1 [0]
      b = (m1 [0] - m0 [0]) / 2.0
      boundaries.append ((m1 [0] + m0 [0]) / 2.0)

      #set first boundary to be before first fromlist point
      if i == 0:
        boundaries [0] = m0 [0] - (boundaries [1] - m0 [0])

    #set last boundary to be after last fromList point
    boundaries.append (fromList [-1][0] + (boundaries [-1] - fromList [-1][0]))

    assert len (boundaries) == len (fromList) + 1

    #return a list of measurements from "toMatches" that are within fromList [fromIndex]'s boundaries
    #the boundaries work like :
    # <= ..... <
    #the search starts at toList [toIndex] and toIndex is incremented globally
    #returning [] means no points found
    def collectMeasurements (fromIndex):
      nonlocal boundaries
      nonlocal toIndex
      nonlocal toList

      measurements = []
      while toIndex < len (toList) and toList [toIndex][0] < boundaries [fromIndex + 1]:
        m = toList [toIndex]
        if m [0] > boundaries [fromIndex]:
          measurements.append (m)
        toIndex += 1
      return measurements

    #return closest (by time) measurement in "measurements" to "m"
    def closest (measurements, measurement):
      if not len (measurements):
        return None
      which = measurements [0]
      distance = 100000000
      for m in measurements [1:]:
        dtemp = math.fabs (measurement [0] - which [0])
        if dtemp < distance:
          which = m
          distance = dtemp
      return which

    #now go down fromList building a pair of lists of measurements matching measurements
    #from fromLists to measurements in toList

    m0 = []
    m1 = []
    toIndex = 0
    for i in range (len (fromList)):
      if toIndex == len (toList):
        break
      measurements = collectMeasurements (i)
      match = closest (measurements, fromList [i])
      if match:
        m0.append (fromList [i])
        m1.append (match)
    return answer (m0, m1)

  #return a tuple of (min, max) of min and max values of measurement values and the min and max
  #values in the minmax argument, if its not None
  def minMaxValue (self, minmax, func):
    filtered = list (filter (func, self.measurements))
    if len (filtered) == 0: return minmax
    small = min ([t[1] for t in filtered])
    big = max ([t[1] for t in filtered])
    if not minmax: return (small, big)
    return (min (small, minmax [0]), max (big, minmax [1]))
  def minMaxUtfSeconds (self, minmax, func):
    filtered = list (filter (func, self.measurements))
    if len (filtered) == 0: return minmax
    small = min ([t[0] for t in filtered])
    big = max ([t[0] for t in filtered])
    if not minmax: return (small, big)
    return (min (small, minmax [0]), max (big, minmax [1]))

#read the next days worth of data from a eso ambient server dump file
def readDay (file, otherDays=None):
  day = Day ()
  s = utils.readNextLine (file)
  if not s: return None
  numMeasurements = int (s.strip())
  assert numMeasurements > 0 and numMeasurements < 1000000
  while numMeasurements > 0:
    fields = file.readline ().strip().split (' ')
    assert len (fields) == 2
    day.addMeasurement (float (fields [0]), float (fields [1]))
    numMeasurements -= 1

  #if otherDays was specified, reset the periodZeroUtfSeconds
  #we have to wait until we've set up day because it needs
  #its julian day to be set
  #periodZeroUtfSeconds is for timestamp filtering based on a start-offset from some other
  #set of measurements
  if otherDays:
    otherDay = otherDays.getDay (day)
    if (otherDay):
      day.periodZeroUtfSeconds = otherDay.measurements [0][0]
      day.otherDayWasFound = True
  else:
    day.otherDayWasFound = True
  return day


#a list of days, read from a file
#otherDays is another set of days to be used to
#normalize period filters with
class Days (object):

  #startDate,stopDate are datetimes
  def __init__(self, path, startDate=None, stopDate=None, otherDays=None):
    self.days = []
    file = open (path, "r")
    day = readDay (file, otherDays)

    #startDate and stopDate are actually datetimes, convert to dates so we
    #can compare just dates
    if startDate and stopDate:
      startDate = startDate.date ()
      stopDate = stopDate.date ()

    while day:
      date = day.datetime.date ()
      if not startDate or (date >= startDate and date <= stopDate):
        self.days.append (day)
      day = readDay (file, otherDays)
    file.close ()

    #make a lookup table mapping dayOfYear to day
    self.dayOfYearToDay = {}
    for day in self.days:
      self.dayOfYearToDay [day.dayOfYear] = day

    #make a lookup table mapping julianDay to day
    self.julianDayToDay = {}
    for day in self.days:
      self.julianDayToDay [day.julianDay] = day

  #return day with the same month/day as that of "day"
  #if our own range of days is in the same calendar year
  #then years are ignored. this lets day matching occur
  #across years. otherwise do an exact date match
  def getDay (self, day):

    if (len (self.days) == 0):
      return None

    if self.days [0].year == self.days [-1].year:
      return self.dayOfYearToDay.get (day.dayOfYear)
    else:
      return self.julianDayToDay.get (day.julianDay)

  def isLeapYear (self):
    if not len (self.days): return False
    return self.days [0].isLeapYear ()

  #set goodness of each day to be that of each day in the Days object "otherDays"
  #days with no corresponding day in "otherDays" are not changed
  def setGoodness (self, otherDays):
    for d in self.days:
      otherDay = otherDays.getDay (d)
      if otherDay:
        d.setGoodness (otherDay)

def matchDays (days, otherDays):
  """
  match days w/ otherDays returning 2 lists, w/ None gaps in the 2nd list if necc.
  :param days: list of Day objects
  :param otherDays: list of Day objects
  :return: tuple of 2 lists, days from each array, matching...
  """
  otherIndex = 0
  days0 = []
  days1 = []
  index = 0

  #if either of the day lists cross a year boundary then we do the matching based on dates
  #otherwise do it based on day of year


  if not len (otherDays) or (days [0].year != days [-1].year or otherDays [0].year != otherDays [-1].year):
    func = lambda d:d.julianDay
  else:
    func = lambda d:d.dayOfYear

  for d in days:
    while otherIndex < len (otherDays) and func (otherDays [otherIndex]) < func (d):
      otherIndex += 1
    if otherIndex == len (otherDays) or func (otherDays [otherIndex]) != func (d):
      days0.append (d)
      days1.append (None)
      continue
    assert func (otherDays [otherIndex]) == func (d)
    days0.append (d)
    days1.append (otherDays [otherIndex])
  return (days0, days1)


def pearsonFromMeasurements (measurements1, measurements2, minPoints):
  """
  compute pearson c.c. between 2 lists of measurements. a sensible result is a number
  between -1 and 1.

  before computing pc, look at the data and test for pathologies and return and out-of-range
  pc value

  2 - lists were empty
  3 - lists had fewer than minPoints items
  4 - exactly 1 list had all identical values, and it was list1
  5 - exactly 1 list had all identical values, and it was list2
  6 - both lists had all identical values
  7 - the lists were identical
  8 - pearson returned Nan for some reason

  :param measurements1:
  :param measurements2:
  :param minPoints
  :return: pc
  """

  assert len (measurements1) == len (measurements2)

  #get just the values from the measurements
  values1 = [m [1] for m in measurements1]
  values2 = [m [1] for m in measurements2]

  if (len (values1) == 0): return 2
  if (len (values1) < minPoints): return 3

  same1 = values1.count (values1[0]) == len (values1)
  same2 = values2.count (values2[0]) == len (values2)

  if same1 and same2: return 6
  if same1: return 4
  if same2: return 5

  same = True
  for v1, v2 in zip (values1, values2):
    if v1 != v2:
      same = False
      break

  if same: return 7

  scipy.seterr (all='raise')
  try:
    pc = scipy.stats.pearsonr (values1, values2) [0]
    return pc
  except FloatingPointError:
    return 8
  finally:
    scipy.seterr (all='warn')



def makeCorrelationMatrix (days, correlations):
  """
  return a 2d matrix of correlations from a list of correlations
  :param days: list of days
  :param correlations: list of (day, day, pc) tuples
  :return: list indexed like list [day-index][day-index] of pc values where the indexes are indexes into "days"
  """
  #make a dictionary of day-to-index mappings so we can do a fast lookup
  indexes = {}
  for index, day in enumerate (days):
    indexes [day] = index

  #make a 2d matrix to fill the pc values in from the correlation list
  m = []
  for i in range (len (days)):
    row = []
    for j in range (len (days)):
      row.append (0)
    m.append (row)

  #fill in the matrix
  for tuple in correlations:
    m [indexes [tuple [0]]][indexes [tuple [1]]] = tuple [2]
    m [indexes [tuple [1]]][indexes [tuple [0]]] = tuple [2]

  return m

def selfCorrelate (days, startSeconds, stopSeconds, pastMidnight, minPoints):
  """
  compute pearson coefficient between every day in "days" with every day in "days"
  :param days: list of Day objects
  :param startSeconds: start of range to accept measurements in seconds
  :param stopSeconds: stop of range to accept measurements in seconds
  :param pastMidnight: True if start/stop are measured past midnight, otherwise they are
   measured past the first measurement timestamp in each of the 2 days
  :return: a list of 4-tuples: (day, day, pc, matchingData), unsorted,
  where matchingData is a 3-tuple: (len of first day's data, len 2nd day's data, len of matching data
  """
  correlations = []
  for twoDays in itertools.combinations (days, 2):
    #get 2 lists of corresponding measurements, 1 for each day in the combination
    f0, f1 = [d.periodFilter (startSeconds, stopSeconds, pastMidnight) for d in twoDays]
    m0, m1 = twoDays [0].matchMeasurements (twoDays [1], f0, f1)
    pc = pearsonFromMeasurements (m0, m1, minPoints)
    correlations.append ((twoDays [0], twoDays[1], pc, (len (twoDays[0].measurements), len (twoDays[1].measurements), len (m0))))
  return correlations

def pairCorrelate (days, otherDays, startSeconds, stopSeconds, pastMidnight, minPoints):
  """
  compute pearson coefficient between every day in "days" with every day in "otherDays"
  these sets of days are assumed to be for the same date range (or close enough)
  :param days: list of Day objects
  :param otherDays: Days object
  :param startSeconds: start of range to accept measurements in seconds
  :param stopSeconds: stop of range to accept measurements in seconds
  :param pastMidnight: True if start/stop are measured past midnight, otherwise they are
   measured past the first measurement timestamp in each of the 2 days
  :return: a list of 4-tuples: (day, day, pc, matchingData), unsorted,
  where matchingData is a 3-tuple: (len of first day's data, len 2nd day's data, len of matching data
  """
  correlations = []
  for twoDays in itertools.combinations (days, 2):
    #get 2 lists of corresponding measurements, 1 for each day in the combination
    day = twoDays [0]
    otherDay = otherDays.getDay (twoDays [1])
    if not otherDay:
      continue
    f0, f1 = [d.periodFilter (startSeconds, stopSeconds, pastMidnight) for d in (day, otherDay)]
    m0, m1 = day.matchMeasurements (otherDay, f0, f1)
    pc = pearsonFromMeasurements (m0, m1, minPoints)
    correlations.append ((day, otherDay, pc, (len (day.measurements), len (otherDay.measurements), len (m0))))
  return correlations


def crossCorrelate (days, otherDays, startSeconds, stopSeconds, pastMidnight, minPoints):
  """
   for each Day object in "days" that has a matching day in "otherDays" compute
   the pc.
  :param days: list of Day objects sorted in ascending julian day
  :param otherDays: another list of Day objects, sorted in ascending julian day
  :param startSeconds: start of range to accept measurements in
  :param stopSeconds: stop of range to accept measurements in
  :param pastMidnight: True if start/stop are measured past midnight, otherwise they are
   measured past the first measurement timestamp in each of the 2 days
  :param minPoints minimum number of points in a set of measurements for pc computation to be valid
  :return: list of 3-tuples: (day, otherDay, pc)
  """
  #get a list of Day objects from otherDays that corresponds to
  #the "days" list. non-matches are None in the resulting list
  days, otherDays = matchDays (days, otherDays)
  correlations = []
  for day, otherDay in zip (days, otherDays):
    if otherDay and day.good and otherDay.good:
      #get 2 lined up lists of measurements suitable for correlation
      f0, f1 = [d.periodFilter (startSeconds, stopSeconds, pastMidnight) for d in (day, otherDay)]
      m0, m1 = day.matchMeasurements (otherDay, f0, f1)
      pc = pearsonFromMeasurements (m0, m1, minPoints)
      correlations.append ((day, otherDay, pc))
    else:
      correlations.append ((day, None, 9))

  return correlations


def goodSquareCorrelate (days, otherDays, minPoints):
  """
   for each Day object in "days" that has a matching day in "otherDays" compute
   the pc. the points used for the pc calculation are within the "good" range
   of each day, and each day's goodness has already been computed.
   days and otherDays have had their empty days padded and are for the same
   date range.
  :param days: list of Day objects sorted in ascending julian day
  :param otherDays: another list of Day objects, sorted in ascending julian day
  :param minPoints minimum number of points in a set of measurements for pc computation to be valid
  :return: 2d array of pc values, array [x][y] is days (x), otherDays.days [y]
  """
  assert len (days) == len (otherDays)
  for day, otherDay in zip (days, otherDays):
    assert day.dayOfYear
    assert day.dayOfYear == otherDay.dayOfYear

  matrix = []
  for i in range (len (days)):
    matrix.append ([])
    for j in range (len (otherDays)):
      matrix [i].append (11)

  for i in range (len (days)):
    day = days [i]
    assert day
    for j in range (len (otherDays)):
      otherDay = otherDays [j]
      assert otherDay
      if not day.good or not otherDay.good:
        continue
      #get 2 lined up lists of measurements suitable for correlation
      f0, f1 = [d.goodPeriodFilter () for d in (day, otherDay)]
      m0, m1 = day.matchMeasurements (otherDay, f0, f1)
      pc = pearsonFromMeasurements (m0, m1, minPoints)
      matrix [i][j] = pc

  return matrix

def goodCorrelate (days, otherDays, minPoints):
  """
   for each Day object in "days" that has a matching day in "otherDays" compute
   the pc. the points used for the pc calculation are within the "good" range
   of each day, and each day's goodness has already been computed
  :param days: list of Day objects sorted in ascending julian day
  :param otherDays: another list of Day objects, sorted in ascending julian day
  :param minPoints minimum number of points in a set of measurements for pc computation to be valid
  :return: list of 3-tuples: (day, otherDay, pc)
  """
  #get a list of Day objects from otherDays that corresponds to
  #the "days" list. non-matches are None in the resulting list
  days, otherDays = matchDays (days, otherDays)
  correlations = []
  for day, otherDay in zip (days, otherDays):
    if otherDay and day.good and otherDay.good:
      #get 2 lined up lists of measurements suitable for correlation
      f0, f1 = [d.goodPeriodFilter () for d in (day, otherDay)]
      m0, m1 = day.matchMeasurements (otherDay, f0, f1)
      pc = pearsonFromMeasurements (m0, m1, minPoints)
      correlations.append ((day, otherDay, pc))
    else:
      correlations.append ((day, None, 9))

  return correlations

def readJulianDayLines (file):
  """
  return a list of strings, one string per line in the file
  where all the strings are on the same julian day as the
  first line read. when this returns
  the file is at the end or the first line of the
  next julian day. and each line in the file starts
  with a julian day
  :param file:
  :return: list of line strings or None if the file had no more lines
  """

  def julianDayFromLine (line):
    fields = line.strip().split (' ')
    assert len (fields) >= 1
    julianDate = float (fields [0])
    fraction = math.modf (julianDate)
    return fraction [1]

  #remember where we are in the file
  pos = file.tell ()

  #get julian day of first line
  line = utils.readNextLine (file)
  if not line: return None
  day = julianDayFromLine (line)

  #go back so we can do the iteration
  file.seek (pos)

  #collect lines into "lines" until the day changes
  lines = []
  while line and julianDayFromLine (line) == day:
    lines.append (line)
    pos = file.tell ()
    line = utils.readNextLine (file)

  #restore file pos to start of new day
  file.seek (pos)
  return lines

def plotSensorSeeingPC (path, correlationsArray, yTicks, seeingDays, site):
  """

  :param path: file to write
  :param correlationsArray: list of list of correlations (day, day, pc)
  :param list of y tick labels
  :param seeingDays: list of Day objects for seeing
  :param site: name of site
  :return:nothing

  """
  assert len (correlationsArray) == len (yTicks)

  m = []
  for correlations in correlationsArray:
    pcs = []
    for day in range (len (seeingDays)):
      if not seeingDays [day]:
        pcs.append (10)
      else:
        assert seeingDays [day] == correlations [day][0]
        pcs.append (correlations [day][2])
    m.append (pcs)

  pdfPages = PdfPages (path)
  fig = plt.figure()
  ax = fig.add_subplot(111)
  cmap = plt.get_cmap ('jet')
  cmap.set_over ('white')
  cs = ax.imshow(m, interpolation='nearest',origin='lower',cmap=cmap, vmax=1.0, extent=[0, len (seeingDays) - 1, 0, len (seeingDays) -1])
  ax.set_title ("PC Seeing & Sensors {:%m/%d/%Y} - {:%m/%d/%Y}\n{}".format (seeingDays [0].datetime, seeingDays [-1].datetime, site))
  ax.set_xlabel ("Day".title ())
  ax.set_ylabel ("Sensor".title ())
  ax.set_yticklabels (yTicks)
  ax.set_yticks (np.arange (0, len (seeingDays) - 1, len (seeingDays) / len (correlationsArray)))
  cbar = fig.colorbar (cs)
  cbar.ax.set_ylabel ("pearson coefficient".title ())
  pdfPages.savefig (fig)
  plt.close (fig)
  pdfPages.close ()


def dayListFromRange (days, startDate, stopDate, leapYear=False):
  """
  return list of Day objects based on taking out a subrange specified by args
  start_date/stop_date take precedence over start_day/stop_day. for subranges
  specified by start/stop date pad the list with Nones where days are missing
  :param days:
  :param startDate: start datetime
  :param stopDate: stop datetime
  :param leapYear: if true insert a padding day at feb 29 if there isn't a day there
  :return: list of Day objects with padding for missing days
  """
  assert startDate and stopDate
  #append toDay to dayList, padding with None's between fromDate up to toDay
  oneDay = datetime.timedelta (days=1)
  def appendBads (dayList, fromDate, uptoDate):
    nonlocal oneDay
    while fromDate < uptoDate:
      if (len (dayList)):
        assert dayList [-1].datetime.date () != fromDate
      badDay = Day (fromDate)
      assert badDay.datetime.date () == fromDate
      dayList.append (badDay)
      fromDate = fromDate + oneDay
  def appendDay (dayList, fromDate, toDay):
    appendBads (dayList, fromDate, toDay.datetime.date ())
    if (len (dayList)):
      assert dayList [-1].datetime.date () != toDay.datetime.date ()
    dayList.append (toDay)

  padded = []
  dayList = [d for d in days.days if d.datetime >= startDate and d.datetime <= stopDate]
  if (len (dayList) == 0):
    appendBads (padded, startDate.date (), stopDate.date () + oneDay)
  else:
    previousDate = None
    for d in dayList:
      if not previousDate:
        appendDay (padded, startDate.date (), d)
      else:
        assert previousDate != d.datetime.date ()
        appendDay (padded, previousDate + oneDay, d)
      previousDate = d.datetime.date ()
    assert previousDate
    appendBads (padded, previousDate + oneDay, stopDate.date () + oneDay)

    #if the day list should have a leap day
    #then we may need to add it

    if leapYear and (not days.isLeapYear ()) and padded [0].dayOfYear < Day.leapDay and padded [-1].dayOfYear > Day.leapDay:
      which = 0
      while padded [which].dayOfYear < Day.leapDay:
        which += 1
      assert (which > 0)
      padded.insert (which, Day (False, dayOfYear=Day.leapDay))

  return padded
