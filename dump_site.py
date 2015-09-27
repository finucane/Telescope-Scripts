#!/usr/local/bin/python3

import os
import ambient
import utils

parser = ambient.ArgumentParser (description="""
print out some info from a site folder specified by --indir
""")

#parse user's command line arguments
args = parser.parse_args ()

def dumpDir (path):
  print ("{}".format (os.path.basename (path)))
  for p in ambient.parameterBasenames:
    days = ambient.Days (os.path.join (path, p + ".orig"), args.start_date, args.stop_date)
    if len (days.days):
      print ("{: >6} {:%m/%d/%y} {:%m/%d/%y} {:>6} days".format (p, days.days[0].datetime, days.days[-1].datetime, len (days.days)))
      for day in days.days:
        if (len (day.measurements) == 0):
          print ("empty day")
        else:
          for m in day.measurements:
            if (len (m) == 0):
              print ("empty measurement")
    else:
      print ("{: >6} empty".format (p))

  print ('\n')

dumpDir (args.indir)

