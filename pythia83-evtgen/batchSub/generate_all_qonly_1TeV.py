import subprocess
import sys

nev = int(sys.argv[1])
kernel = int(sys.argv[2])
CF = float(sys.argv[3])
CA = float(sys.argv[4])

prefs = ['1000']
#prefs = ['0100']
exps = ['0000','0011','1000','1011']

evt_per = 10000
njobs = nev // evt_per

config_dir = "/global/homes/s/sambt/Jets/optimal-classifiers/pythia83-evtgen/DIRE-config-files/"
base_cfg_quark = "-n "+str(evt_per)+" -k "+str(kernel)+" --cf "+str(CF)+" --ca "+str(CA) + " -c "+config_dir+"H2qqbar-ee-m1000.cmnd"+" -c "+config_dir+"fixGluonRecoil.cmnd"
base_cfg_gluon = "-n "+str(evt_per)+" -k "+str(kernel)+" --cf "+str(CF)+" --ca "+str(CA) + " -c "+config_dir+"H2gg-ee-m1000.cmnd"+" -c "+config_dir+"fixGluonRecoil.cmnd"

""""
#Generating baseline H --> gg samples
for p in prefs:
    name = 'H2gg-CF{0:.1f}CA{1:.1f}-sqg{2}'.format(CF,CA,p)
    cfg_str = base_cfg_gluon + " -c "+config_dir+"splitting-sqg{0}.cmnd".format(p)
    for i in range(1,njobs+1):
        command = ["sbatch","run_job.sh",name,str(i),str(kernel),cfg_str]
        subprocess.run(command)

#Generating H --> qq samples with exponential splitting variants
for p in prefs:
    for ex in exps:
        if ex == '0000':
            name = 'H2qq-CF{0:.1f}CA{1:.1f}-sqg{2}'.format(CF,CA,p)
            cfg_str = base_cfg_quark + " -c "+config_dir+"splitting-sqg{0}.cmnd".format(p)
        else:
            name = 'H2qq-CF{0:.1f}CA{1:.1f}-sqg{2}-esq{3}'.format(CF,CA,p,ex)
            cfg_str = base_cfg_quark + " -c " + config_dir + "splitting-sqg{0}-esq{1}.cmnd".format(p,ex)
        for i in range(1,njobs+1):
            command = ["sbatch","run_job.sh",name,str(i),str(kernel),cfg_str]
            subprocess.run(command)
"""
#Generating H --> qq samples with no g --> gg splittings
for p in prefs:
    for ex in exps:
        name = 'H2qq-m1000-CF{0:.1f}CA{1:.1f}-sq{2}-esq{3}'.format(CF,CA,p,ex)
        cfg_str = base_cfg_quark + " -c " + config_dir + "q_only/splitting-sqg{0}-esq{1}.cmnd".format(p,ex)
        for i in range(1,njobs+1):
            command = ["sbatch","run_job.sh",name,str(i),str(kernel),cfg_str]
            subprocess.run(command)
