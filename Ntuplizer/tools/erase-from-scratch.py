import os
from optparse import OptionParser


parser = OptionParser()
parser.add_option('--wn', action="store",type="int",dest="WN",default=0)
(options, args) = parser.parse_args()

wn = []

if not options.WN == 0:
  print "Erasing all files on t3wn%i" %options.WN
  wn.append(options.WN)
else:  
  print "Erasing all files on all t3 wns"
  for i in range(10,41):
    wn.append(i)
    
print ""  
for w in wn:
 cmd =   "qrsh -q debug.q -l hostname=t3wn%i ls /scratch/$USER" %w
 print ""
 print "Listing all files in your scratch on t3wn%i:"%w
 os.system(cmd)
 cmd = "qrsh -q debug.q -l hostname=t3wn%i find /scratch/$USER/ -user $USER -exec rm -rf {} \\;" %w
 print ""
 print "Removing all files in t3wn%i /scratch/$USER"%w
 os.system(cmd)