import subprocess

import sys

name = str(sys.argv[1])
nev = int(sys.argv[2])
kernel = int(sys.argv[3])
configs = sys.argv[4:]

config_dir = "/global/homes/s/sambt/Jets/optimal-classifiers/pythia83-evtgen/DIRE-config-files/"

CFs = [0.0001,1/3,2/3,1,4/3,5/3,2,7/3,8/3,3]
CAs = [3.0]
evt_per = 10000
njobs = nev // evt_per

print("---------- Running CF/CA Scan ----------")
print("****************************************")
for CF in CFs:
    for CA in CAs:
        print("Submitting jobs for {0} events with CF = {1}, CA = {2}, configs = {3}".format(nev,CF,CA,configs))

        cfg_str = "-n "+str(evt_per)+" --cf "+str(CF)+" --ca "+str(CA)+" -k "+str(kernel)
        for c in configs:
            cfg_str += " -c "+config_dir+c+".cmnd"
        print("config string : "+cfg_str)
            
        for i in range(1,njobs+1):
            command = ["sbatch","-C","haswell","-q","regular","-t","10","run_job.sh",name+"-CF{0:.1f}CA{1:.1f}".format(CF,CA),str(i),str(kernel),cfg_str]
            subprocess.run(command)
        print("****************************************")
            
            
