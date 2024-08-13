import uproot as ur
filename=f"lambda_recon.root"
events = ur.open(f'{filename}:events')
arrays = events.arrays()

#get the truth value of theta*
px=arrays_sim["MCParticles.momentum.x"][:,2]
py=arrays_sim["MCParticles.momentum.y"][:,2]
pz=arrays_sim["MCParticles.momentum.z"][:,2]
tilt=-0.025
pt=np.hypot(px*np.cos(tilt)-pz*np.sin(tilt), py)
theta=np.arctan2(pt,pz*np.cos(tilt)+px*np.sin(tilt))

#compute the value of theta* using the clusters in the ZDC
xc=arrays["HcalFarForwardZDCClusters.position.x"]
yc=arrays["HcalFarForwardZDCClusters.position.y"]
zc=arrays["HcalFarForwardZDCClusters.position.z"]
E=arrays_sim["HcalFarForwardZDCClusters.energy"]

rc=np.sqrt(xc**2+yc**2+zc**2)
xcp=xc*np.cos(tilt)-zc*np.sin(tilt)
ycp=yc
zcp=zc*np.cos(tilt)+xc*np.sin(tilt)

E=arrays_sim["HcalFarForwardZDCClusters.energy"][i]

px_recon,py_recon,pz_recon=np.sum(E*xcp/rc, axis=-1),np.sum(E*ycp/rc, axis=-1),np.sum(E*zcp/rc, axis=-1)
pt_recon=np.hypot(px_recon,py_recon)
theta_recon=np.arctan2(pt_recon, np.sum(E, axis=-1))



