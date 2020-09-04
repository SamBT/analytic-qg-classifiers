import subprocess
import sys

jobname = str(sys.argv[1])
filename = str(sys.argv[2])
process = str(sys.argv[3])
nev = int(sys.argv[4])
kernel = int(sys.argv[5])

#Cvals = [0.0001,0.3333,0.6666,1.0,1.3333,1.6666,2.0,2.3333,2.6666,3.0]
#Cvals = [0.3333]
CFs = [0.0001,1/3,2/3,1,4/3,5/3,2,7/3,8/3,3]
CAs = [3.0]
evt_per = 10000
njobs = nev // evt_per
remain = nev % evt_per

print("---------- Running CF/CA Scan ----------")
print("****************************************")
for CF in CFs:
    for CA in CAs:
        print("Submitting jobs for {0} events with CF = {1}, CA = {2}, process = {3}".format(nev,CF,CA,process))
        for i in range(1,njobs+1):
            command = ["sbatch","-C","haswell","-q","regular","run_job.sh",jobname+"_CF_{0:.1f}_CA_{1:.1f}".format(CF,CA),filename+"_CF_{0:.1f}_CA_{1:.1f}".format(CF,CA),str(i),process,str(evt_per),str(CF),str(CA),str(kernel)]
            subprocess.run(command)
        #run remainder of jobs, if needed
        if (remain > 0):
            command = ["sbatch","-C","haswell","-q","regular","run_job.sh",jobname+"_CF_{0:.1f}_CA_{1:.1f}".format(CF,CA),filename+"_CF_{0:.1f}_CA_{1:.1f}".format(CF,CA),str(i),process,str(remain),str(CF),str(CA),str(kernel)]
            subprocess.run(command)
        print("****************************************")
            
            
