import subprocess
import datetime

minenergy = 10 # in GeV
maxenergy = 100 # in GeV
#pid = "neutron"
#name = "neutron_{}_{}GeV".format(minenergy,maxenergy)
#pid = "proton"
#name = "proton_{}_{}GeV".format(minenergy,maxenergy)
#pid = "kaon0L"
#name = "kaon0_{}_{}GeV".format(minenergy,maxenergy)
#pid = "kaon+"
#name = "kaon_{}_{}GeV".format(minenergy,maxenergy)
#pid = "pi+"
#name = "pi_{}_{}GeV".format(minenergy,maxenergy)
#pid = "gamma"
#name = "gamma_{}_{}GeV".format(minenergy,maxenergy)
#pid = "e-"
#name = "ele_{}_{}GeV".format(minenergy,maxenergy)
pid = "pi0"
name = "pi0_{}_{}GeV".format(minenergy,maxenergy)
batch = 10
queue = 20
begin = 0
num_cpu=2
mem=5.6
##kcms05 might be problem
num_punch=0
with open("punch.txt","r") as f:
  num_punch=len(f.readlines())
script="""Universe = vanilla
getenv = True
should_transfer_files=yes
when_to_transfer_output= ON_EXIT
request_memory = {mem}GB
request_cpus = {num_cpu}
transfer_input_files = build
Executable = rungps.sh
Requirements = ( machine !="kcms05" )
transfer_output_files = box
output = condor/{name}_$(Process).out
error = condor/{name}_$(Process).error
Log = condor/{name}_$(Process).log
Arguments = run_gps $(Process) {begin} box/{name} {pid} {minenergy} {maxenergy} {batch} {num_cpu}
Queue {queue}""".format(name=name,pid=pid,minenergy=minenergy,maxenergy=maxenergy,batch=batch,begin=begin,queue=queue,num_cpu=num_cpu,mem=mem)

f=open("buf_gps.co","w")
f.write(script)
f.close()

import subprocess
print("#script below")
print(script)
a=raw_input("Do you want run this? y/n : ")
if(a=="y"):
  subprocess.call(["condor_submit","--batch-name","{}_{}_{}".format(name,queue,num_punch),"buf_gps.co"])
  print("job submitted")
  with open("punch.txt","a") as f:
    f.write("{}_{}_{};{}\n".format(name,queue,num_punch,str(datetime.datetime.now())))
