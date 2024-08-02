import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from mpl_toolkits.mplot3d import Axes3D
import uproot
import pandas as pd
import awkward as ak
import matplotlib.colors as mcolors
import math
#from pylorentz import Momentum4
import scipy
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
import os
import vector

plt.figure()
hep.set_style(hep.style.CMS)
hep.set_style("CMS")
plt.close()

def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))


tree = uproot.open(f'sim_output/semi_coherent/{os.environ["DETECTOR_CONFIG"]}_semi_coherent.eicrecon.tree.edm4eic.root')['events']
cluster_nhits = tree['EcalFarForwardZDCClusters/EcalFarForwardZDCClusters.nhits'].array()
cluster_energy = tree['EcalFarForwardZDCClusters/EcalFarForwardZDCClusters.energy'].array()
cluster_theta = tree['EcalFarForwardZDCClusters/EcalFarForwardZDCClusters.intrinsicTheta'].array()
cluster_phi = tree['EcalFarForwardZDCClusters/EcalFarForwardZDCClusters.intrinsicPhi'].array()
cluster_x = tree['EcalFarForwardZDCClusters/EcalFarForwardZDCClusters.position.x'].array()  
cluster_y = tree['EcalFarForwardZDCClusters/EcalFarForwardZDCClusters.position.y'].array()
cluster_z = tree['EcalFarForwardZDCClusters/EcalFarForwardZDCClusters.position.z'].array()

cluster_nhits = pd.DataFrame({'cluster_nhits': np.array(cluster_nhits.tolist(), dtype=object)})
cluster_energy = pd.DataFrame({'cluster_energy': np.array(cluster_energy.tolist(), dtype=object)})
cluster_theta = pd.DataFrame({'cluster_theta': np.array(cluster_theta.tolist(), dtype=object)})
cluster_phi = pd.DataFrame({'cluster_phi': np.array(cluster_phi.tolist(), dtype=object)})
cluster_x = pd.DataFrame({'cluster_x': np.array(cluster_x.tolist(), dtype=object)})
cluster_y = pd.DataFrame({'cluster_y': np.array(cluster_y.tolist(), dtype=object)})
cluster_z = pd.DataFrame({'cluster_z': np.array(cluster_z.tolist(), dtype=object)})

reco_energy = tree['EcalFarForwardZDCRecHits/EcalFarForwardZDCRecHits.energy'].array()
reco_x = tree['EcalFarForwardZDCRecHits/EcalFarForwardZDCRecHits.position.x'].array()  
reco_y = tree['EcalFarForwardZDCRecHits/EcalFarForwardZDCRecHits.position.y'].array()
reco_z = tree['EcalFarForwardZDCRecHits/EcalFarForwardZDCRecHits.position.z'].array()

reco_energy = pd.DataFrame({'reco_energy': np.array(reco_energy.tolist(), dtype=object)})
reco_x = pd.DataFrame({'reco_x': np.array(reco_x.tolist(), dtype=object)})
reco_y = pd.DataFrame({'reco_y': np.array(reco_y.tolist(), dtype=object)})
reco_z = pd.DataFrame({'reco_z': np.array(reco_z.tolist(), dtype=object)})

mc_id = tree['MCParticles/MCParticles.PDG'].array()
mc_status = tree['MCParticles/MCParticles.generatorStatus'].array()
mc_x = tree['MCParticles/MCParticles.momentum.x'].array()  
mc_y = tree['MCParticles/MCParticles.momentum.y'].array()
mc_z = tree['MCParticles/MCParticles.momentum.z'].array()

mc_id = pd.DataFrame({'mc_id': np.array(mc_id.tolist(), dtype=object)})
mc_status = pd.DataFrame({'mc_status': np.array(mc_status.tolist(), dtype=object)})
mc_x = pd.DataFrame({'mc_x': np.array(mc_x.tolist(), dtype=object)})
mc_y = pd.DataFrame({'mc_y': np.array(mc_y.tolist(), dtype=object)})
mc_z = pd.DataFrame({'mc_z': np.array(mc_z.tolist(), dtype=object)})

dg = pd.concat([cluster_nhits, cluster_energy, cluster_theta, cluster_phi, cluster_x, cluster_y, cluster_z,reco_energy,reco_x,reco_y,reco_z,mc_id,mc_status,mc_x,mc_y,mc_z],axis=1)

def rotateY(xdata, zdata, angle):
    s = np.sin(angle)
    c = np.cos(angle)
    rotatedz = c*zdata - s*xdata
    rotatedx = s*zdata + c*xdata
    return rotatedx, rotatedz

Mp = 0.938 #GeV
mom = 110 #GeV/c/nucleon

z = 92
a = 238

cluster_energy = np.concatenate(dg.cluster_energy)
X = np.concatenate(dg.cluster_x)
Y = np.concatenate(dg.cluster_y)
Z = np.concatenate(dg.cluster_z)
X, Z = rotateY(X,Z,0.025)
cluster_theta = np.arccos(Z/np.sqrt(Z**2+Y**2+X**2))
cluster_phi = np.arctan2(Y,X)
cluster_eta = -np.log(np.tan(cluster_theta/2))

#cluster_vec = Momentum4.e_m_eta_phi(cluster_energy,0,cluster_eta,cluster_phi)*1000
#cluster_vec_boost = cluster_vec.boost_particle(Momentum4(np.sqrt(Mp**2+mom**2) ,0,0,-mom))

pz = cluster_energy * np.cos(cluster_theta)
pt = cluster_energy * np.sin(cluster_theta)
px = pt * np.cos(cluster_phi)
py = pt * np.sin(cluster_phi)
cluster_theta = np.arccos(Z/np.sqrt((X**2+Y**2+Z**2)))*1000

cluster_vec = vector.array({"px": px, "py": py, "pz": pz, "E": cluster_energy})*1000
boost = vector.array({"px": [0], "py": [0], "pz": [-mom], "energy": [np.sqrt(Mp**2+mom**2)]})
cluster_vec_boost = cluster_vec.boost_p4(boost)

X = np.concatenate(dg.mc_x)
Y = np.concatenate(dg.mc_y)
Z = np.concatenate(dg.mc_z)
id = np.concatenate(dg.mc_id)

X, Z = rotateY(X,Z,0.025)
mc_energy = np.sqrt(np.array(X)**2 + np.array(Y)**2 + np.array(Z)**2)
mc_theta = np.arccos(Z/np.sqrt(Z**2+Y**2+X**2))
mc_phi = np.arctan(Y/X)
mc_eta = -np.log(np.tan(mc_theta/2))
mc_theta = np.arccos(Z/np.sqrt((X**2+Y**2+Z**2)))*1000

#mc_vec = Momentum4.e_m_eta_phi(mc_energy,0,mc_eta,mc_phi)*1000
#mc_vec_boost = mc_vec.boost_particle(Momentum4(np.sqrt(Mp**2+mom**2) ,0,0,-mom))

mc_vec = vector.array({"px": X, "py": Y, "pz": Z, "energy": mc_energy})*1000
mc_vec_boost = mc_vec.boost_p4(boost)

fig, ax = plt.subplots(2,1,sharex=True,figsize=(15,15))
plt.sca(ax[0])
#plt.hist2d(cluster_vec_boost[0],cluster_theta,bins=(np.logspace(-2,1,20),np.linspace(0,5,100)),label=f'Reconstructed',cmin=0.01,cmap="YlOrBr")
plt.scatter(mc_vec_boost.energy[id==22],mc_theta[id==22],c='blue')
plt.scatter(cluster_vec_boost.energy,cluster_theta,c='orange')
plt.title(f'Fragment: A = {a}, Z = {z}\n Boost Reference Frame: U238 Beam Rest Frame\nReconstructed at ZDC Ecal')
#plt.colorbar()
'''
Energy = []
Theta = []
Phi = []
for i in E.index:
    x = Momentum4(np.array(pE[i])-1,np.array(ppx[i])-1,np.array(ppy[i])-1,np.array(ppz[i])-1)*1000
    x = x.boost_particle(Momentum4(np.sqrt(Mp**2+mom**2) ,0,0,-mom))
    j = np.arccos((np.array(ppz[i])-1)/np.sqrt((np.array(ppx[i])-1)**2+(np.array(ppy[i])-1)**2+(np.array(ppz[i])-1)**2))*1000
    plt.scatter(x[0],j,c='b')
    Energy.append(x[0])
    Theta.append(j)
    Phi.append(np.tan(np.sqrt((np.array(ppy[i])-1)/np.sqrt((np.array(ppx[i])-1)))))
plt.scatter(x[0],j,c='b',label='BeAGLE')
'''
plt.ylabel('Theta (mRad)')
plt.xscale('log')
plt.ylim(0,4)
#plt.xlim(1E-5,1E-2)


plt.sca(ax[1])
hist1, bins1 = np.histogram(cluster_vec_boost.energy[cluster_theta<4],bins=np.logspace(-2,1,100))
bins1 = bins1[1:]/2 + bins1[:-1]/2
plt.errorbar(bins1,hist1/sum(hist1),yerr=np.sqrt(hist1)/sum(hist1),label=f'Reconstructed Clusters',fmt='-o',c='orange')

condition1 = id==22
condition2 = mc_theta<4
hist, bins = np.histogram(mc_vec_boost.energy[condition1 & condition2],bins=np.logspace(-2,1,100))
bins = bins[1:]/2 + bins[:-1]/2
plt.errorbar(bins,hist/sum(hist),yerr=np.sqrt(hist)/sum(hist),label=f'MC Particles',fmt='-o',c='b')

plt.xscale('log')
#plt.yscale('log')
plt.ylabel('Count (Normalized)')
plt.xlabel('Energy (MeV)')

plt.legend()
plt.xlim(1E-2,1E1)
fig.tight_layout(pad=-0.5)

with PdfPages(f'results/{os.environ["DETECTOR_CONFIG"]}/semi_coherent/plots.pdf') as pdf:
    pdf.savefig(fig)
