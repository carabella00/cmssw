import os
import sys
from itertools import islice
submit = 'universe = vanilla\n' ##writing .sub file
submit += 'arguments = "$(argument)"\n'
submit += 'output = submit01.out\n'
submit += 'error = submit01.err\n'
submit += 'log = submit01.log\n'
submit += '+JobFlavour = "tomorrow"\n' ##finish writing .sh file
submit += 'queue\n'
submitName = 'submit01.sub'
sub1 = open(submitName,'w')
sub1.write(submit+'\n')
sub1.close()
nfile = 10 ##number of files per job
filename = '../../TrackFindingTMTT/test/MCsamples/1040/RelVal/TTbar/PU200.txt' ## you need grid certificate
with open(filename,'r') as f:
    counter1 = 1 ## the nth job
    while True:
        lines = list(islice(f, nfile))
        if not lines:
            break
        counter2 = 1 ## the nth file in one job
        for line in lines:
            if counter2 == 1:
               input='root://cms-xrd-global.cern.ch//'+line.rstrip() ## yon don't need 'root:..' if you don't have grid certificate
            else:
               input =input+',root://cms-xrd-global.cern.ch//'+line.rstrip()
            counter2+=1
        create = '#!/bin/bash\n' ##writng .sh file
        create += 'export CMSSW_PROJECT_SRC=/afs/cern.ch/user/j/jingyan/CMS_TEST/CMSSW_10_6_0/src\n'
        create += 'cd $CMSSW_PROJECT_SRC\n'
        create += 'eval `scramv1 runtime -sh`\n'
        create += 'export X509_USER_PROXY=/afs/cern.ch/user/j/jingyan/x509up_u122075\n' ## you don't need this line if you don't have grid certificate
        create += 'cd /afs/cern.ch/user/j/jingyan/CMS_TEST/CMSSW_10_6_0/src/L1Trigger/TrackFindingTracklet/test\n'
        create += 'cmsRun L1TrackNtupleMakerTestParsing_cfg.py inputFiles='+input+' outputFile='+str(100*nfile)+'events_D21_Hybrid'+str(counter1)+'.root  maxEvents='+str(100*nfile)+'\n'
        createName = 'submit'+str(counter1)+'.sh'
        sub2 = open(createName,'w')
        sub2.write(create+'\n')
        sub2.close()
        counter1+=1
        os.system('chmod 755 '+createName) ## make .sh file executable
        os.system('condor_submit '+ submitName+' executable='+createName) ## submit the job






