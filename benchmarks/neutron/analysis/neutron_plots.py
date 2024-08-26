import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl, awkward as ak, sys, uproot as ur
import mplhep as hep
hep.style.use("CMS")

plt.rcParams['figure.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'
plt.rcParams['savefig.bbox']='tight'

plt.rcParams["figure.figsize"] = (7, 7)
config=sys.argv[1].split("/")[1]  #results/{config}/neutron
outdir=sys.argv[1]+"/"
try:
    import os
    os.mkdir(outdir[:-1])
except:
    pass

#read files
arrays_sim={}
for p in 20,30,40,50,60,70,80:
    arrays_sim[p] = ur.open(f'sim_output/neutron/{config}_rec_neutron_{p}GeV.edm4hep.root:events')\
                    .arrays()

#get reconstructed theta, eta, and E
def gauss(x, A,mu, sigma):
    return A * np.exp(-(x-mu)**2/(2*sigma**2))

#get the truth pseudorapidity and energy
for array in arrays_sim.values():
    tilt=-0.025
    px=array['MCParticles.momentum.x'][:,2]
    py=array['MCParticles.momentum.y'][:,2]
    pz=array['MCParticles.momentum.z'][:,2]
    p=np.sqrt(px**2+py**2+pz**2)
    
    pxp=px*np.cos(tilt)-pz*np.sin(tilt)
    pyp=py
    pzp=pz*np.cos(tilt)+px*np.sin(tilt)
    
    array['E_truth']=np.hypot(p, 0.9406)
    array['eta_truth']=1/2*np.log((p+pzp)/(p-pzp))
    array['theta_truth']=np.arccos(pzp/p)

#weighted avg of theta of cluster centers, by energy
for array in arrays_sim.values():
    tilt=-0.025
    x=array['HcalEndcapPInsertClusters.position.x']
    y=array['HcalEndcapPInsertClusters.position.y']
    z=array['HcalEndcapPInsertClusters.position.z']
    E=array['HcalEndcapPInsertClusters.energy']
    r=np.sqrt(x**2+y**2+z**2)
    
    xp=x*np.cos(tilt)-z*np.sin(tilt)
    yp=y
    zp=z*np.cos(tilt)+x*np.sin(tilt)
    
    w=E
    
    array['theta_recon']=np.sum(np.arccos(zp/r)*w, axis=-1)/np.sum(w, axis=-1)
    array['eta_recon']=-np.log(np.tan(array['theta_recon']/2))
    
    
    array['E_Hcal']=np.sum(array['HcalEndcapPInsertClusters.energy'], axis=-1)#*20/12.5
    array['E_Ecal']=np.sum(array['EcalEndcapPInsertClusters.energy'], axis=-1)#*27/20
    
    array['E_recon']=array['E_Hcal']/(0.70-.30/np.sqrt(array['E_Hcal']))\
                        +array['E_Ecal']/(0.82-0.35/np.sqrt(array['E_Ecal']))

#plot theta residuals:
from scipy.optimize import curve_fit

fig, axs=plt.subplots(1,2, figsize=(16,8))
plt.sca(axs[0])
p=40
eta_min=3.4; eta_max=3.6
y,x,_=plt.hist(1000*(arrays_sim[p]['theta_recon']-arrays_sim[p]['theta_truth'])\
               [(arrays_sim[p]['eta_truth']>eta_min)&(arrays_sim[p]['eta_truth']<eta_max)], bins=50,
                    range=(-10,10), histtype='step')
bc=(x[1:]+x[:-1])/2
slc=abs(bc)<3
# try:
fnc=gauss
sigma=np.sqrt(y[slc])+(y[slc]==0)
p0=(100, 0, 5)
coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,sigma=list(sigma))
xx=np.linspace(-5,5,100)
plt.plot(xx,fnc(xx,*coeff))
# except:
#     pass
plt.xlabel("$\\theta_{rec}-\\theta_{truth}$ [mrad]")
plt.ylabel("events")
plt.title(f"$p={p}$ GeV, ${eta_min}<\\eta<{eta_max}$")

r=[3.0,3.2,3.4,3.6,3.8]
for eta_min, eta_max in zip(r[:-1],r[1:]):
    xvals=[]
    sigmas=[]
    dsigmas=[]
    for p in 20,30,40, 50, 60, 70, 80:
        y,x=np.histogram(1000*(arrays_sim[p]['theta_recon']-arrays_sim[p]['theta_truth'])\
                         [(arrays_sim[p]['eta_truth']>eta_min)&(arrays_sim[p]['eta_truth']<eta_max)],
                         bins=50, range=(-10,10))
        bc=(x[1:]+x[:-1])/2
        slc=abs(bc)<3
        fnc=gauss
        p0=(100, 0, 5)
        #print(bc[slc],y[slc])
        sigma=np.sqrt(y[slc])+(y[slc]==0)
        try:
            coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,sigma=list(sigma))
            sigmas.append(np.abs(coeff[2]))
            dsigmas.append(np.sqrt(var_matrix[2][2]))
            xvals.append(p)
        except:
            pass
    plt.sca(axs[1])
    plt.errorbar(xvals, sigmas, dsigmas, ls='', marker='o', label=f"${eta_min}<\\eta<{eta_max}$")
plt.xlabel("$p_{n}$ [GeV]")
plt.ylabel("$\\sigma[\\theta]$ [mrad]")
plt.ylim(0)
plt.legend()
plt.tight_layout()
plt.savefig(outdir+"neutron_theta_recon.pdf")

fig, axs=plt.subplots(1,2, figsize=(16,8))
plt.sca(axs[0])
p=50
eta_min=3.4; eta_max=3.6
y,x,_=plt.hist(((arrays_sim[p]['E_recon']-arrays_sim[p]['E_truth'])/arrays_sim[p]['E_truth'])\
               [(arrays_sim[p]['eta_truth']>eta_min)&(arrays_sim[p]['eta_truth']<eta_max)&(arrays_sim[p]['E_Hcal']>0)], bins=50,
                    range=(-.5,.5), histtype='step')
bc=(x[1:]+x[:-1])/2
slc=abs(bc)<0.4
# try:
p0=(100, 0, 0.3)
fnc=gauss
sigma=np.sqrt(y[slc])+(y[slc]==0)
coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,sigma=list(sigma))
xx=np.linspace(-.5,.5,100)
plt.plot(xx,fnc(xx,*coeff))
# except:
#     pass
plt.xlabel("$(E_{rec}-E_{truth})/E_{truth}$")
plt.ylabel("events")
plt.title(f"$p={p}$ GeV, ${eta_min}<\\eta<{eta_max}$")

r=[3.0,3.2,3.4,3.6,3.8]
for eta_min, eta_max in zip(r[:-1],r[1:]):
    xvals=[]
    sigmas=[]
    dsigmas=[]
    for p in 20,30,40, 50, 60, 70, 80:
        y,x=np.histogram(((arrays_sim[p]['E_recon']-arrays_sim[p]['E_truth'])/arrays_sim[p]['E_truth'])\
                       [(arrays_sim[p]['eta_truth']>eta_min)&(arrays_sim[p]['eta_truth']<eta_max)], bins=50,  range=(-.5,.5))
        bc=(x[1:]+x[:-1])/2
        slc=abs(bc)<0.5
        fnc=gauss
        p0=(100, 0, 0.3)
        #print(bc[slc],y[slc])
        sigma=np.sqrt(y[slc])+(y[slc]==0)
        try:
            coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,sigma=list(sigma))
            sigmas.append(np.abs(coeff[2]))
            dsigmas.append(np.sqrt(var_matrix[2][2]))
            xvals.append(p)
        except:
            pass
    plt.sca(axs[1])
    plt.errorbar(xvals, sigmas, dsigmas, ls='', marker='o', label=f"${eta_min}<\\eta<{eta_max}$")
plt.xlabel("$p_{n}$ [GeV]")
plt.ylabel("$\\sigma[E]/E$")
plt.ylim(0)
plt.legend()
plt.tight_layout()
plt.savefig(outdir+"neutron_energy_recon.pdf")
