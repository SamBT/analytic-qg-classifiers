import subprocess
import sys

process = str(sys.argv[1])
nev = int(sys.argv[2])
kernel = int(sys.argv[3])

CFs = [0.0001,1/3,2/3,1,4/3,5/3,2,7/3,8/3,3]
CAs = [3.0]

print("---------- Running CF/CA Scan ----------")
print("****************************************")
for CF in CFs:
    for CA in CAs:
        print("Running {0} events with CF = {1}, CA = {2}, process = {3}".format(nev,CF,CA,process))
        command = ["./GenEvts",process+".cmnd",str(nev),str(CF),str(CA),str(kernel)]
        subprocess.run(command)
        command = ["mv","test.root","output-files/"+process+"/kernel{0}/".format(kernel)+process+"-CF-{0:.2f}-CA-{1:.2f}".format(CF,CA)+".root"]
        subprocess.run(command)
        print("****************************************")
