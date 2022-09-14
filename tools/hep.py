import argparse
import subprocess
import datetime,sys
parser = argparse.ArgumentParser()
parser.add_argument('-tw','--towerwidth',type=int, default=2, help='towerwidth*2 + 1 tower range particle gun incident')
parser.add_argument('-pn','--particlename',type=str, default="", help='particle name for particle gun')
parser.add_argument('-e','--energy',type=str, default="70", help='energy for particle gun')
parser.add_argument('-s','--seed',type=int, default=-1, help='random seed')
#python3 hep.py -pn e- -e 50 -s 123 -tw 2
args=parser.parse_args()
max_energy=-1
if(args.particlename!=""):
  pid=args.particlename
  if("-" in args.energy):
      energy=int(args.energy.split("-")[0])
      max_energy=int(args.energy.split("-")[1])
  else:
    energy=int(args.energy)
  seed=args.seed
else:
  print(["uu","dd","cc","ss","gg"])
  pid=raw_input("Pick one of these particles : ")
  seed=raw_input("Set seed number : ")
jets=[1,2,3,4,5,6,21]
pdgid={"dd":1,"uu":2,"ss":3,"cc":4,"bb":5,"tt":6,"gg":21,
        "e-":11,"gamma":22,"pi0":111,"pi+":211,"kaon0L":130,"kaon+":321}
P8 = "P8generic"
mod="" #mod*2 +1 tower width will cover
isjet=False
batch = 10
queue = 100
begin = 0
num_cpu=1
mem=5.3
if(pid in pdgid.keys()):
  P8 = "P8ptcgen"
  cmnd = "qqgun"
  if(max_energy==-1):name = "{}_{}GeV".format(pid,energy)
  else:name = "{}_{}-{}GeV".format(pid,energy,max_energy)
  if(abs(pdgid[pid]) in jets):
      isjet=True
  else:
      batch = 20
      mod=str(args.towerwidth)
else:
  if(pid=="wwqqep"):
    cmnd = "wwqqep"
    name = "ee_to_WW_250GeV"
  if(pid=="zzqqep"):
    cmnd = "zzqqep"
    name = "ee_to_ZZ_250GeV"

if(cmnd=="qqgun"):
  with open('../build/Gen/qqgun.cmnd','r') as f:
    c=f.readlines()
    with open('../build/Gen/{}.cmnd'.format(name),'w') as d:
      for line in c:
        if("nEvt_HERE" in line):
          line=line.replace("nEvt_HERE",str(batch))
        if("PID_HERE" in line):
          line=line.replace("PID_HERE",str(pdgid[pid]))
        if("E_HERE" in line):
          line=line.replace("E_HERE",str(energy))
        if(isjet):
          if("pT_HAT_MIN_HERE" in line):
            line=line.replace("pT_HAT_MIN_HERE",str(30))
          if("COLOR_SINGLET_HERE" in line):
            line=line.replace("COLOR_SINGLET_HERE","false")
        else:
          if("pT_HAT_MIN_HERE" in line):
            line=line.replace("pT_HAT_MIN_HERE",str(0))
          if("COLOR_SINGLET_HERE" in line):
            line=line.replace("COLOR_SINGLET_HERE","true")
        d.write(line)
  #cmnd="{}".format(name)

with open('../build/DRsim/run_hepmc.mac','r') as f:
  c=f.readlines()
  with open('../build/DRsim/run_{}.mac'.format(name),'w') as d:
    for line in c:
      if("nEvt_HERE" in line):
        line=line.replace("nEvt_HERE",str(batch))
      d.write(line)

with open('runhep.sh','r') as f:
  c=f.readlines()
  with open('run_buf.sh','w') as d:
    for line in c:
      d.write(line)
punch=[]
with open("punch.txt","r") as f:
  punch_cards=f.readlines()
  for line in punch_cards:
    if("#"==line[0]):
      continue
    punch.append(line)

is_on=0
ex_seed=[]
if(seed==-1):
  for i in range(len(punch_cards)):
      if("#"==punch_cards[i][0]):
        buf_seed=int(punch_cards[i].split(";")[0].split("_")[-1])
        if(is_on!=-1 and name in punch_cards[i]):
          ex_seed.append(buf_seed)
      else:
        if(name in punch_cards[i]):
          is_on=-1
          break
  print(is_on)
  if(is_on==-1):
    sys.exit(1)
  for i in range(max(ex_seed)):
    if(not i in ex_seed):
      is_on=i
      break
  if(is_on!=0):
      seed=is_on
  else:
      seed=max(ex_seed)+1

begin = queue*seed
##occationally few job stuck in loop
#Requirements = ( machine !="kcms23" )
script="""Universe = vanilla
getenv = True
should_transfer_files=yes
when_to_transfer_output= ON_EXIT
request_memory = {mem}GB
request_cpus = {num_cpu}
transfer_input_files = ../build
Executable = run_buf.sh
transfer_output_files = box
output = condor/{name}_$(Process).out
error = condor/{name}_$(Process).error
Log = condor/{name}_$(Process).log
Arguments = {name} $(Process) {begin} box/{name} {batch} {P8} {max_energy} {mod}
Queue {queue}""".format(name=name,batch=batch,begin=begin,queue=queue,num_cpu=num_cpu,mem=mem,P8=P8,max_energy=max_energy,mod=mod)

f=open("buf_hep.co","w")
f.write(script)
f.close()

import subprocess
print("#script below")
print(script)
#a=raw_input("Do you want run this? y/n : ")
a="y"
if(a=="y"):
  subprocess.call(["condor_submit","--batch-name","{}_{}_{}".format(name,queue,seed),"buf_hep.co"])
  print("job submitted")
  with open("punch.txt","a") as f:
    f.write("{}_{}_{};{}\n".format(name,queue,seed,str(datetime.datetime.now())))
