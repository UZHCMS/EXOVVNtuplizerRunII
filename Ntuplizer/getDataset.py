import subprocess, os, sys
from optparse import OptionParser, OptionValueError

prefix = 'dcap://t3se01.psi.ch:22125/'


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


cdir = os.getcwd()


import datetime

d_today = datetime.datetime.now()

d_today = str(d_today).split('.')[0]
d_today = d_today.replace(' ', '-').replace(':','')

print d_today



usage = "usage: python getDataset.py"
parser = OptionParser(usage)

#/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/

parser.add_option("-p", "--path", default="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/", type="string", help="path", dest="path")
parser.add_option("-c", "--chunk", default=1, type="int", help="chunk", dest="chunk")
parser.add_option("-a", "--analysis", default="BsTauTau", type="string", help="analysis channel", dest="analysis")
parser.add_option("-t", "--type", default="mc", type="string", help="type", dest="type")

(options, args) = parser.parse_args()

print 'Path = ', options.path
print 'Chunks = ', options.chunk


jobdir = cdir + '/job/' + options.analysis + '_' + options.type + '_' + d_today

ensureDir(jobdir)


#RJpsi_ParkingBPH1_2020-01-19-130139_20200119130138_BsPhiTauTau_mutau/ParkingBPH1/ParkingBPH1_Run2018A-05May2019-v1_20200119130138_BsPhiTauTau_mutau/200119_120740/0000

out = subprocess.Popen(['ls', options.path], 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.STDOUT)


listoffiles = []



for line in iter(out.stdout.readline,''):

   line = line.rstrip()
 
#   print line 

   listoffiles.append(prefix + '/' + options.path + '/' + line)
  
#   if len( line.split())!=9: continue

#   ndir = line.rstrip().split()[8]

#   print 'Adding ...', prefix + '/' + options.path + '/' + ndir

#   out2 = subprocess.Popen(['uberftp', '-ls', prefix + '/' + options.path + '/' + ndir],
#                           stdout=subprocess.PIPE, 
#                           stderr=subprocess.STDOUT)



#sys.exit(1)


for file in listoffiles:
    print file


print 
print len(listoffiles), 'files are detected'

#sys.exit(1)

#listoffiles = list(chunks(listoffiles,options.chunk))

print len(listoffiles), 'chunks are created'


for ijob, filename in enumerate(listoffiles):
    print ijob, filename

#    if ijob==2: break
#    input = ','.join(files)

    # Now, write the shell script into job directory
    
    jobscript = jobdir + '/job_' + str(ijob) + '.sh'
    outfile = jobdir + '/flatTuple_' + str(ijob) + '.root'

    os.system("cp job_template.sh " + jobscript)

    
    with open(jobscript) as f:
        data_lines = f.read()
        
    data_lines = data_lines.replace('INFILE', filename).replace('OUTFILE', outfile)
        
    with open(jobscript, mode="w") as f:
        f.write(data_lines)


    command = 'sbatch -p wn --account=t3 --output=' + jobdir + '/out.' + str(ijob) + ' ' + jobscript
    os.system(command)
