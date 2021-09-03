import numpy as np
import ROOT as rt
import datetime
from tqdm import tqdm
from particle import PDGID

now=datetime.datetime.now()
_fb=56*3 # raw image pixel number
#_cb=8
_cb=0 # image padding
num_point=2048
name_file=["ele.root","pi.root"] #image root file path
name_cls=["el","pi"] #output file name
save_path=["npzs"] #output directory path
for nf,nc in zip(name_file,name_cls):
  event=rt.TChain("event")
  event.Add(nf)
  images=[]
  points=[]
  ptcs=[]
  pt_gen=[]
  E_Gen=[]
  width_Gen=[]
  SC_ratio=[]
  dr_ecorr=[]
  print(event.GetEntries())
  for i in tqdm(range(int(event.GetEntries()))):
      event.GetEntry(i)
      pt_gen.append(event.pt_Gen)
      E_Gen.append(event.E_Gen)
      width_Gen.append(event.width_Gen)
      SC_ratio.append(event.E_S/event.E_C)
      dr_ecorr.append(event.E_DRcorr)
      images.append([np.array(list(event.image_ecor_s)[:_fb*_fb]).reshape((_fb,_fb))[_cb:_fb-_cb,_cb:_fb-_cb],
          np.array(list(event.image_ecor_c)[:_fb*_fb]).reshape((_fb,_fb))[_cb:_fb-_cb,_cb:_fb-_cb],
          np.array(list(event.image_n_s)[:_fb*_fb]).reshape((_fb,_fb))[_cb:_fb-_cb,_cb:_fb-_cb],
          np.array(list(event.image_n_c)[:_fb*_fb]).reshape((_fb,_fb))[_cb:_fb-_cb,_cb:_fb-_cb]])#s c ns nc
      point=np.array([np.array(list(event.fiber_phi)),
          np.array(list(event.fiber_theta)),
          np.array(list(event.fiber_depth))/1000.,
          np.array(list(event.fiber_ecor_s)),
          np.array(list(event.fiber_ecor_c)),
          ]).transpose()#phi theta depth s c
      gen=np.array([np.array(list(event.ptc_phi)),
          np.array(list(event.ptc_theta)),
          np.array(list(event.ptc_p)),
          np.array(list(event.ptc_E)),
          np.array(list(event.ptc_pid)),
          np.zeros(len(event.ptc_phi)),
          np.zeros(len(event.ptc_phi)),
          np.zeros(len(event.ptc_phi)),
          np.zeros(len(event.ptc_phi)),
          ]).transpose()#phi,theta,p,E,px,py,pz,ch,nh,cl,nl
          #np.array(list(event.ptc_px)),
          #np.array(list(event.ptc_py)),
          #np.array(list(event.ptc_pz)),
      point=sorted(point, key=lambda pnt:pnt[3],reverse=True)
      #print(len(point),len(gen))
      #print([len(i) for i in gen])
      if(len(point)<2048):
          point=np.concatenate([point,np.zeros((2048-len(point),5))],axis=0)
      if(len(point)>2048):
          point=point[:2048]
      points.append(point)

      gen=sorted(gen, key=lambda pnt:pnt[3],reverse=True)
      if(len(gen)<64):
          gen=np.concatenate([gen,np.zeros((64-len(gen),9))],axis=0)
      if(len(gen)>64):
          gen=gen[:64]
      for j in range(len(gen)):
        is_l=-1
        if(PDGID(gen[j][4]).is_hadron==0):is_l=0
        if(PDGID(gen[j][4]).is_lepton==0):is_l=1
        if(is_l!=-1):gen[j][5+2*is_l+int(PDGID(gen[j][4]).charge==0)]=1

      ptcs.append(gen)


  print("depth",np.mean(list(event.fiber_depth)))
  images=np.array(images,dtype='float32')
  #print([len(i) for i in point])
  points=np.array(points,dtype='float32')
  ptcs=np.array(ptcs,dtype='float32')
  print(len(images),len(point))
  #print(images.shape,point.shape)
  np.savez("{}/dr_{}".format(save_path,nc),image=images,point=points,ptc=ptcs,pt_gen=pt_gen,dr_ecorr=dr_ecorr,E_Gen=E_Gen,width_Gen=width_Gen,SC_ratio=SC_ratio)

  del event,images,points
  print(nf,nc,"time costs ",datetime.datetime.now()-now)

print("time costs ",datetime.datetime.now()-now)
