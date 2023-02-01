import os
import numpy as np
import uproot as ur
import awkward as ak
import matplotlib.pyplot as plt
import matplotlib as mpl
import mplhep 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--rec_file', type=str, help='Reconstructed track file.')
parser.add_argument('--config', type=str, help='Momentum configuration.')
parser.add_argument('--nevents', type=float, help='Number of events to process.')
parser.add_argument('--results_path', type=str, help='Output directory.')
args = parser.parse_args()
kwargs = vars(args)

rec_file = args.rec_file
Nevents = int(args.nevents)
r_path = args.results_path + '/truth_reconstruction/' #Path for output figures and file.
Dconfig = 'epic' + args.config.split('_epic')[1].strip() #Detector config
config = args.config.split('_epic')[0].strip()

for array in ur.iterate(rec_file + ':events',['MCParticles/MCParticles.generatorStatus',
                                        'MCParticles/MCParticles.PDG',
                                        'MCParticles/MCParticles.momentum.x',
                                        'MCParticles/MCParticles.momentum.y',
                                        'MCParticles/MCParticles.momentum.z',
                                        'ReconstructedChargedParticles/ReconstructedChargedParticles.PDG',
                                        'ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.x',
                                        'ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.y',
                                        'ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.z',
                                        'ReconstructedChargedParticlesAssociations/ReconstructedChargedParticlesAssociations.simID',
                                        'ReconstructedChargedParticlesAssociations/ReconstructedChargedParticlesAssociations.recID'],step_size=Nevents):
    PDG_mc = array['MCParticles/MCParticles.PDG']  #Monte Carlo (MC) particle numbering scheme.
    px_mc = array['MCParticles/MCParticles.momentum.x']
    py_mc = array['MCParticles/MCParticles.momentum.y']
    pz_mc = array['MCParticles/MCParticles.momentum.z']
    PDG_rc = array['ReconstructedChargedParticles/ReconstructedChargedParticles.PDG']
    px_rc = array['ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.x']
    py_rc = array['ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.y']
    pz_rc = array['ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.z']
    simID = array['ReconstructedChargedParticlesAssociations/ReconstructedChargedParticlesAssociations.simID']
    recID = array['ReconstructedChargedParticlesAssociations/ReconstructedChargedParticlesAssociations.recID']
    #SimID and recID contain the indices of the MCParticles and ReconstructedParticles entry for that event.

### MCParticles Variables
momentum_mc = np.sqrt(((px_mc**2)+(py_mc**2)+(pz_mc**2)))
theta_mc = np.arctan2(np.sqrt(px_mc**2+py_mc**2), pz_mc)
phi_mc = np.arctan2(py_mc, px_mc)
### ReconstructedParticles Variables
momentum_rc = np.sqrt(((px_rc**2)+(py_rc**2)+(pz_rc**2)))
theta_rc = np.arctan2(np.sqrt(px_rc**2+py_rc**2), pz_rc)
phi_rc = np.arctan2(py_rc, px_rc)

booll = (PDG_mc[simID])==(PDG_rc[recID]) #boolean that allows events where the same particle is reconstructed
boolean_pion = np.logical_or(ak.Array(PDG_mc[simID][booll])==-211, ak.Array(PDG_mc[simID][booll])==+211) #boolean that allows events involving pions
boolean_proton = np.logical_or(ak.Array(PDG_mc[simID][booll])==-2212, ak.Array(PDG_mc[simID][booll])==+2212) #boolean that allows events involving protons
boolean_electron = ak.Array(PDG_mc[simID][booll])==11 #boolean that allows events involving electrons
boolean_neutron = ak.Array(PDG_mc[simID][booll])==2112 #boolean that allows events involving neutrons
boolean_photon = ak.Array(PDG_mc[simID][booll])==22 #boolean that allows events involving photons

### MCParticles variables list
MC_list = [ak.Array(momentum_mc[simID][booll]),                     #Momentum
           ak.Array(theta_mc[simID][booll]),                        #Theta
           ak.Array(phi_mc[simID][booll]),                          #Phi
           -np.log(np.tan((ak.Array(theta_mc[simID][booll]))/2))]   #Eta
### ReconstructedParticles variables list
RC_list = [ak.Array(momentum_rc[recID][booll]),                     #Momentum
           ak.Array(theta_rc[recID][booll]),                        #Theta
           ak.Array(phi_rc[recID][booll]),                          #Phi
           -np.log(np.tan((ak.Array(theta_rc[recID][booll]))/2))]   #Eta
title_list = ['Momentum','Theta','Phi','Eta']
### MC Momentum for different particles list
M_list = [ak.Array(momentum_mc[simID][booll]),
          ak.Array(momentum_mc[simID][booll][boolean_pion]),
          ak.Array(momentum_mc[simID][booll][boolean_proton]),
          ak.Array(momentum_mc[simID][booll][boolean_electron]),
          ak.Array(momentum_mc[simID][booll][boolean_neutron]),
          ak.Array(momentum_mc[simID][booll][boolean_photon])]

#Marker Size in plots
if Nevents == 100:
    ssize = 1
else:
    ssize = 0.01

#Particle type for Single events
particle = config.split('-')[0].strip()
particle_dict = {'e':[boolean_electron,'Electrons'],'pi':[boolean_pion,'Pions']}


####################################################################################################
    #Ratio 
####################################################################################################

for i in range(len(MC_list)): #Repeat the following steps for each variable (momentum,theta,phi,eta)
    MCparts = MC_list[i] #MCParticles events to be plotted on x-axis
    RCparts = RC_list[i] #ReconstructedParticles events
    X_list = [ak.Array(MCparts),
          ak.Array(MCparts[boolean_pion]),
          ak.Array(MCparts[boolean_proton]),
          ak.Array(MCparts[boolean_electron]),
          ak.Array(MCparts[boolean_neutron]),
          ak.Array(MCparts[boolean_photon])]
    Y_list = [ak.Array(RCparts),
          ak.Array(RCparts[boolean_pion]),
          ak.Array(RCparts[boolean_proton]),
          ak.Array(RCparts[boolean_electron]),
          ak.Array(RCparts[boolean_neutron]),
          ak.Array(RCparts[boolean_photon])]
    X_plot = list(np.zeros(len(X_list)))
    Y_plot = list(np.zeros(len(X_list)))

    for j in range(len(X_list)): #Repeat the following steps for each particle (pions,protons,electrons,neutrons,photons)
        X = X_list[j]
        Y = Y_list[j]
        X_len = ak.count(X,axis=None)
        Y_len = ak.count(Y,axis=None)
        if X_len > Y_len: 
            F_boolean = np.ones_like(Y) == 1
        else: 
            F_boolean = np.ones_like(X) == 1
        X_s = np.array(ak.flatten(X[F_boolean])) #Filtered lists
        Y_s = np.array(ak.flatten(Y[F_boolean]))
        if i == 0:   #Momentum
            ratio = np.array((ak.Array(Y_s)/ak.Array(X_s)))
        else: #Angle difference
            ratio = np.array((ak.Array(Y_s)-(ak.Array(X_s))))
        X_plot[j] = X_s
        Y_plot[j] = ratio

    fig = plt.figure()
    gs = fig.add_gridspec(3, 2, wspace=0)
    (ax1, ax2), (ax3, ax4),(ax5, ax6) = gs.subplots(sharex=True, sharey=True)
    # fig.suptitle('')
    if i == 1:   # for theta
        X_plot[0],X_plot[1],X_plot[2],X_plot[3],X_plot[4],X_plot[5] = -X_plot[0],-X_plot[1],-X_plot[2],-X_plot[3],-X_plot[4],-X_plot[5]
    ax1.scatter(X_plot[0], Y_plot[0], s = ssize)
    ax2.scatter(X_plot[1], Y_plot[1], s = ssize)
    ax3.scatter(X_plot[2], Y_plot[2], s = ssize)
    ax4.scatter(X_plot[3], Y_plot[3], s = ssize)
    ax5.scatter(X_plot[4], Y_plot[4], s = ssize)
    ax6.scatter(X_plot[5], Y_plot[5], s = ssize)

    if i == 0:  # for momentum
        ax1.set_ylabel('rc/mc') #ratio
        ax3.set_ylabel('rc/mc')
        ax5.set_ylabel('rc/mc') 
        title ='ratio'
        ax1.set_yscale('log')
        ax1.set_xscale('log')
    else:       # for angles
        ax1.set_ylabel('rc-mc') #difference
        ax3.set_ylabel('rc-mc')
        ax5.set_ylabel('rc-mc') 
        title ='difference'
    ax2.set_title('Pions')
    ax3.set_title('Protons')
    ax4.set_title('Electrons')
    ax5.set_title('Neutrons')
    ax6.set_title('Photons')

    if i == 3: #Eta
        ax1.set_xlim(-5.5,5.5)
    if i == 1: #Theta
        ax1.set_xlim(-np.pi,0)
        ax5.set_xlabel('- %s mc'%(title_list[i]))
        ax6.set_xlabel('- %s mc'%(title_list[i]))
        x_range = [0,np.pi]
    else:
        ax5.set_xlabel('%s mc'%(title_list[i]))
        ax6.set_xlabel('%s mc'%(title_list[i]))
        x_range = list(ax1.get_xlim())
    fig.set_figwidth(20)
    fig.set_figheight(10)
    ax1.set_title('%s %s  %s  %s events\n DETECTOR_CONFIG: %s'%(title_list[i],title,config,Nevents,Dconfig))
    plt.savefig(os.path.join(r_path, '%s_%s_%s.png' %  (title_list[i],title,config)))
    plt.close()


###################################################################################################
     #Ratio vs momentum
###################################################################################################

    if i > 0: #for each variable theta, phi, and eta
        for j in range(len(M_list)): #Repeat the following steps for each particle (pions,protons,electrons,neutrons,photons)
            X = X_list[j]
            Y = Y_list[j]
            M_mc = M_list[j]
            boolean_M = np.ones_like(M_mc) == 1
            X_s = np.array(ak.flatten(X[boolean_M])) 
            Y_s = np.array(ak.flatten(Y[boolean_M])) 
            M_s = np.array(ak.flatten(M_mc))
            ratio = np.array((ak.Array(Y_s)-(ak.Array(X_s))))
            X_plot[j] = M_s
            Y_plot[j] = ratio

        fig = plt.figure()
        gs = fig.add_gridspec(3, 2, wspace=0)
        (ax1, ax2), (ax3, ax4),(ax5, ax6) = gs.subplots(sharex=True, sharey=True)
        ax1.scatter(X_plot[0], Y_plot[0], s = ssize)
        ax2.scatter(X_plot[1], Y_plot[1], s = ssize)
        ax3.scatter(X_plot[2], Y_plot[2], s = ssize)
        ax4.scatter(X_plot[3], Y_plot[3], s = ssize)
        ax5.scatter(X_plot[4], Y_plot[4], s = ssize)
        ax6.scatter(X_plot[5], Y_plot[5], s = ssize)
        ax1.set_xscale('log')
        ax1.set_ylabel('rc-mc')
        ax3.set_ylabel('rc-mc')
        ax5.set_ylabel('rc-mc')
        ax5.set_xlabel('Momentum mc')
        ax6.set_xlabel('Momentum mc')
        ax2.set_title('Pions')
        ax3.set_title('Protons')
        ax4.set_title('Electrons')
        ax5.set_title('Neutrons')
        ax6.set_title('Photons')
        fig.set_figwidth(20)
        fig.set_figheight(10)
        ax1.set_title('%s Difference Vs Momentum  %s  %s events\n DETECTOR_CONFIG: %s'%(title_list[i],config,Nevents,Dconfig))
        plt.savefig(os.path.join(r_path, '%s_difference_vs_momentum_%s.png' %  (title_list[i],config)))


###################################################################################################
     #Correlation
###################################################################################################

    #Repeat the following steps for each variable (momentum,theta,phi,eta)
    X_len = ak.count(MCparts,axis=None)
    Y_len = ak.count(RCparts,axis=None)
    if X_len > Y_len: 
        F_boolean = np.ones_like(RCparts) == 1
    else: 
        F_boolean = np.ones_like(MCparts) == 1
    X_s = np.array(ak.flatten(MCparts[F_boolean])) 
    Y_s = np.array(ak.flatten(RCparts[F_boolean])) 
 
    #Histogram
    if i == 0 and particle in particle_dict.keys(): #Momentum in Single events
        h, xedges, yedges, image = plt.hist2d(x=X_s,y= Y_s,  bins = 11) 
    else:    
        h, xedges, yedges, image = plt.hist2d(x=X_s,y= Y_s,  bins = 11, range = [x_range,x_range]) 
    plt.close()

    col_sum = ak.sum(h,axis=-1) #number of events in each (verticle) column 
    norm_h = [] #norm_h is the normalized matrix
    norm_h_text = []
    for j in range(len(col_sum)):
        if col_sum[j] != 0:
            norm_c = h[j]/col_sum[j] #normalized column = column values divide by sum of the column
        else:
            norm_c = h[j]
        norm_h.append(norm_c)
        norm_c_text = [ '%.3f' % elem for elem in norm_c ] #display value to 3 dp
        norm_h_text.append(norm_c_text)
    
    fig, axs = plt.subplots(1, 2, figsize=(20, 10))
    if i == 0 and particle in particle_dict.keys(): #Momentum in Single events
        axs[0].hist2d(x=X_s,y=Y_s, bins = 11)
    else: 
        axs[0].hist2d(x=X_s,y=Y_s, bins = 11,range = [x_range,x_range]) 
    mplhep.hist2dplot(H=norm_h,norm=mpl.colors.LogNorm(vmin= 1e-4, vmax= 1),labels=norm_h_text, xbins = xedges, ybins = yedges, ax=axs[1])
    
    axs[0].set_title('%s Histogram'%(title_list[i]))
    axs[0].set_xlabel('%s_mc'%(title_list[i]))
    axs[0].set_ylabel('%s_rc'%(title_list[i]))
    axs[1].set_xlabel('%s_mc'%(title_list[i]))
    axs[1].set_ylabel('%s_rc'%(title_list[i]))
    axs[1].set_title('%s Correlation'%(title_list[i]))
    fig.suptitle('%s  %s events\n DETECTOR_CONFIG: %s'%(config,Nevents,Dconfig))
    plt.savefig(os.path.join(r_path, '%s_correlation_%s.png' %  (title_list[i],config)))


###################################################################################################
     #Phi vs Theta plots
###################################################################################################

def particle_plots(boolean_particle):
    #filtered lists w.r.t the particle
    theta_mc_fil = ak.Array(theta_mc[simID][booll])[boolean_particle]
    theta_rc_fil = ak.Array(theta_rc[recID][booll])[boolean_particle]
    phi_mc_fil = ak.Array(phi_mc[simID][booll])[boolean_particle]
    phi_rc_fil = ak.Array(phi_rc[recID][booll])[boolean_particle]

    theta_mc_fil_len = ak.count(theta_mc_fil,axis=None)
    theta_rc_fil_len = ak.count(theta_rc_fil,axis=None)
    if theta_mc_fil_len > theta_rc_fil_len:
        F_boolean = np.ones_like(theta_rc_fil) == 1
    else: 
        F_boolean = np.ones_like(theta_mc_fil) == 1
    #filtered lists w.r.t length
    theta_mc_F = np.array(ak.flatten(theta_mc_fil[F_boolean]))
    theta_rc_F = np.array(ak.flatten(theta_rc_fil[F_boolean]))
    phi_mc_F = np.array(ak.flatten(phi_mc_fil[F_boolean]))
    phi_rc_F = np.array(ak.flatten(phi_rc_fil[F_boolean]))
    ratio = np.array((ak.Array(theta_rc_F)-(ak.Array(theta_mc_F))))

    fig = plt.figure()
    gs = fig.add_gridspec(2, 2, wspace=0.01)
    (ax1, ax2), (ax3, ax4) = gs.subplots(sharex=True, sharey=True)
    ax1.scatter(-theta_mc_F, ratio, s = ssize)
    ax2.scatter(-theta_rc_F, ratio, s = ssize)
    ax3.scatter(-theta_mc_F, phi_mc_F, s = ssize)
    ax4.scatter(-theta_rc_F, phi_rc_F, s = ssize)
    ax1.set_ylabel('rc-mc')
    ax2.set_ylabel('rc-mc')
    ax3.set_ylabel('Phi mc')
    ax4.set_ylabel('Phi rc')
    ax1.set_xlabel('- Theta mc')
    ax2.set_xlabel('- Theta rc')
    ax3.set_xlabel('- Theta mc')
    ax4.set_xlabel('- Theta rc')
    fig.set_figwidth(20)
    fig.set_figheight(10)

if particle in particle_dict.keys():
    boolean_particle = particle_dict[particle][0]
    particle_name = particle_dict[particle][1]
    particle_plots(boolean_particle)

    plt.suptitle('%s in %s  %s events\n DETECTOR_CONFIG: %s'%(particle_name,config,Nevents,Dconfig))
    plt.savefig(os.path.join(r_path, '%s_%s.png' %  (particle_name,config)))
else:
    for i in [[boolean_photon,'Photons'],[boolean_electron,'Electrons'],[boolean_pion,'Pions']]:
        boolean_particle = i[0]
        particle_name = i[1]
        particle_plots(boolean_particle)

        plt.suptitle('%s in %s  %s events\n DETECTOR_CONFIG: %s'%(particle_name,config,Nevents,Dconfig))
        plt.savefig(os.path.join(r_path, '%s_%s.png' %  (particle_name,config)))





