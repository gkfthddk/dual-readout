import numpy as np
import ROOT as rt
import datetime
now=datetime.datetime.now()
_fb=56*3
#_cb=8
_cb=0
num_point=2048
#name_file=["~/dream/elenshower.root","~/dream/pienshower.root"]
name_file=["~/dream/el.root","~/dream/pi.root","~/dream/ga.root","~/dream/pi0.root"]
name_cls=["el","pi","ga","pi0"]
for nf,nc in zip(name_file,name_cls):
  event=rt.TChain("event")
  event.Add(nf)
  images=[]
  points=[]
  pt_gen=[]
  dr_ecorr=[]
  print(event.GetEntries())
  for i in range(int(event.GetEntries())):
      event.GetEntry(i)
      pt_gen.append(event.pt_Gen/1000.)
      dr_ecorr.append(event.E_DRcorr)
      images.append([np.array(list(event.image_ecor_s)[:_fb*_fb]).reshape((_fb,_fb))[_cb:_fb-_cb,_cb:_fb-_cb],
          np.array(list(event.image_ecor_c)[:_fb*_fb]).reshape((_fb,_fb))[_cb:_fb-_cb,_cb:_fb-_cb],
          np.array(list(event.image_n_s)[:_fb*_fb]).reshape((_fb,_fb))[_cb:_fb-_cb,_cb:_fb-_cb],
          np.array(list(event.image_n_c)[:_fb*_fb]).reshape((_fb,_fb))[_cb:_fb-_cb,_cb:_fb-_cb]])
      point=np.array([np.array(list(event.fiber_phi)),
          np.array(list(event.fiber_theta)),
          np.array(list(event.fiber_depth))/1000.,
          np.array(list(event.fiber_ecor_s)),
          np.array(list(event.fiber_ecor_c)),
          ]).transpose()
      point=sorted(point, key=lambda pnt:pnt[3],reverse=True)
      if(len(point)<2048):
          point=np.concatenate([point,np.zeros((2048-len(point),5))],axis=0)
      if(len(point)>2048):
          point=point[:2048]
      points.append(point)
      #point=[np.array(list(event.fiber_phi)),
      #    np.array(list(event.fiber_theta)),
      #    np.array(list(event.fiber_depth))/1000.,
      #    np.array(list(event.fiber_ecor_s)),
      #    np.array(list(event.fiber_ecor_c)),
      #    ]
      #point=sorted(point, key=lambda pnt:pnt[3],reverse=True)
      #if(len(point[0])<2048):
      #    point=np.concatenate([point,np.zeros((5,2048-len(point[0])))],axis=1)
      #if(len(point[0])>2048):
      #    point=np.array(point)[:,:2048]
      #points.append(np.swapaxes(point,0,1))


  print("depth",np.mean(list(event.fiber_depth)))
  images=np.array(images,dtype='float32')
  print([len(i) for i in point])
  points=np.array(points,dtype='float32')
  print(images.shape,point.shape)
  np.savez("npzs/dr20_{}".format(nc),image=images,point=points,pt_gen=pt_gen,dr_ecorr=dr_ecorr)# save files for each in name_cls
  del event,images,points

print("time costs ",datetime.datetime.now()-now)
