import subprocess
import sys

nev = int(sys.argv[1])
kernel = int(sys.argv[2])
CF = float(sys.argv[3])
CA = float(sys.argv[4])
radius = 0.4
if len(sys.argv) == 6:
    radius = float(sys.argv[5])

exp_z = [0,1]
exp_z2 = [0,1]

if nev < 10000:
    evt_per = nev
    njobs = 1
else:
    evt_per = 10000
    njobs = nev // evt_per

config_dir = "/global/homes/s/sambt/Jets/optimal-classifiers/pythia83-evtgen/DIRE-config-files/"
base_cfg_quark = "-n "+str(evt_per)+" -k "+str(kernel)+" -r "+str(radius)+" --cf "+str(CF)+" --ca "+str(CA) + \
                " -c "+config_dir+"H2qqbar-ee.cmnd"+" -c "+config_dir+"fixGluonRecoil.cmnd"+" -c "+config_dir+"splitting-sqg0100.cmnd"
base_cfg_gluon = "-n "+str(evt_per)+" -k "+str(kernel)+" -r "+str(radius)+" --cf "+str(CF)+" --ca "+str(CA) + \
                " -c "+config_dir+"H2gg-ee.cmnd"+" -c "+config_dir+"fixGluonRecoil.cmnd"+" -c "+config_dir+"splitting-sqg0100.cmnd"

#Generating vanilla H --> gg with the given radius and CF, CA settings
for i in range(1,njobs+1):
    name = "H2gg-CF{0:.1f}CA{1:.1f}-sqg0100-r{2:.1f}".format(CF,CA,radius)
    command = ["sbatch","run_job.sh",name,str(i),str(kernel),base_cfg_gluon]
    subprocess.run(command)

#Generating H --> qq samples with exponential splitting variants
for ez in exp_z:
    for ez2 in exp_z2:
        print("z = {0}, z2 = {1}".format(ez,ez2))
        esq = "00{0}{1}".format(ez,ez2)
        cmd1 = "DireGeneralizedKernel:softExps:Dire_fsr_qcd_1->1&21=0.,0.,{0}.,{1}.".format(ez,ez2)
        cmd2 = "DireGeneralizedKernel:softExps:Dire_fsr_qcd_1->21&1=0.,0.,{0}.,{1}.".format(ez,ez2)
        name = 'H2qq-CF{0:.1f}CA{1:.1f}-sqg0100-esq{2}-r{3:.1f}'.format(CF,CA,esq,radius)
        cfg_str = base_cfg_quark + " -e " + cmd1 + " -e "+cmd2
        for i in range(1,njobs+1):
            command = ["sbatch","run_job.sh",name,str(i),str(kernel),cfg_str]
            subprocess.run(command)
