import subprocess

minenergy = 10 # in GeV
maxenergy = 100 # in GeV
#pid = "neutron"
#name = "neutron_{}GeV".format(energy)
#pid = "kaon0L"
#name = "kaon0_{}GeV".format(energy)
#pid = "pi+"
#pid = "gamma"
#name = "gamma_{}_{}GeV".format(minenergy,maxenergy)
pid = "e-"
name = "ele_{}_{}GeV".format(minenergy,maxenergy)
#pid = "pi0"
#name = "pi0_{}_{}GeV".format(minenergy,maxenergy)
batch = 20
queue = 600
begin = 400
num_cpu=6
mem=6

script="""Universe = vanilla
getenv = True
should_transfer_files=yes
when_to_transfer_output= ON_EXIT
request_memory = {mem}GB
request_cpus = {num_cpu}
Requirements = ( machine !="kcms14" && machine !="kcms13" || machine == "kcms10")
transfer_input_files = build
Executable = rungps.sh
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
  subprocess.call(["condor_submit","--batch-name","{}_{}".format(name,batch),"buf_gps.co"])
  print("job submitted")
