##=============================================================================================
##
##    Script to run Pythia + Dire event generation in an automated way
##
##          Eric M. Metodiev, MIT, 2018 - 2019
##
##    python3 running.py
##
##=============================================================================================

import subprocess

# number of events and center of mass energy (GeV)
nev, ecm = 100000, 200
Cs = [0.0001,0.3333,0.6666,1.0000,1.3333,1.6666,2.0000,2.3333,2.6666,3.0000]

# iterate over process: qq (0) and gg (1)
for process in [0,1]:
  for CA in Cs:
	  for CF in Cs:
	    # command line argument to run event generation
	    cmd = ["./events", "dire.cmnd", "-nev", str(nev), "-proc", str(process), "-ecm", str(ecm), "-CF", str(CF), "-CA", str(CA)]
	    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	    
	    # print output
	    [print(line) for line in p.stdout]
	    p.wait()
