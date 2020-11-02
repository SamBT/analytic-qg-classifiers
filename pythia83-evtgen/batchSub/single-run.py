import subprocess
import sys

name = str(sys.argv[1])
nev = int(sys.argv[2])
kernel = int(sys.argv[3])
CF = float(sys.argv[4])
CA = float(sys.argv[5])
configs = sys.argv[6:]

config_dir = "/global/homes/s/sambt/Jets/optimal-classifiers/pythia83-evtgen/DIRE-config-files/"

evt_per = 10000
njobs = nev // evt_per

print("Submitting jobs for {0} events with CF = {1}, CA = {2}, configs = {3}".format(nev,CF,CA,configs))

cfg_str = "-n "+str(evt_per)+" -k "+str(kernel)+" --cf "+str(CF)+" --ca "+str(CA)
for c in configs:
    cfg_str += " -c "+config_dir+c+".cmnd"
print("config string : "+cfg_str)
for i in range(1,njobs+1):
    command = ["sbatch","-C","haswell","-q","regular","-t","15","run_job.sh",name,str(i),str(kernel),cfg_str]
    subprocess.run(command)

            
            
