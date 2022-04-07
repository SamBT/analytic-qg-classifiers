import subprocess
import sys

nev = int(sys.argv[1])
kernel = int(sys.argv[2])
radius = 0.4
if len(sys.argv) == 6:
    radius = float(sys.argv[5])

CFs = [0.0000001,1/3,2/3,1.0,4/3,5/3,2.0,7/3,8/3,3.0]
CA = 3.0

if nev < 10000:
    evt_per = nev
    njobs = 1
else:
    evt_per = 10000
    njobs = nev // evt_per

config_dir = "/global/homes/s/sambt/Jets/optimal-classifiers/pythia83-evtgen/DIRE-config-files/"

#Generating H --> qq samples with exponential splitting variants
for CF in CFs:
    cfg_str = "-n "+str(evt_per)+" -k "+str(kernel)+" -r "+str(radius)+" --cf "+str(CF)+" --ca "+str(CA) + \
              " -c "+config_dir+"H2qqbar-ee.cmnd"+" -c "+config_dir+"fixGluonRecoil.cmnd"+" -c "+config_dir+"splitting-sqg0100.cmnd"
    name = 'H2qq-CF{0:.1f}CA{1:.1f}-sqg0100-r{2:.1f}'.format(CF,CA,radius)
    for i in range(1,njobs+1):
        command = ["sbatch","run_job.sh",name,str(i),str(kernel),cfg_str]
        subprocess.run(command)
