# Instructions 

This branch you are reading - "diffractive-phi-benchmarks" under physics_benchmark, is currently being developed by Kong.Tu (kongtu@bnl.gov). 

Here the instructions are ONLY intended for people who wants to give it a shot at full ATHENA detector simulation, especially if you do not have much experience in ATHENA (yet). In addition, if you want to follow this particular setup/framework, you need to have ALL of the following:

- BNL RACF account, if you do not have and would like to have one, follow this link, https://wiki.bnl.gov/eic/index.php/Getting_Started. When you just get a new account, you need to setup your EIC environment. If you are desperate enough to use my setup (my setup changes all the time due to other tasks), see /gpfs02/eic/ztu/.cshrc. But I strongly suggest you contact Kolja Kauder <kkauder@gmail.com> (don't tell him you got it from here:) ) 


- Your own ready-to-use HEPMC MC sample.

- Register for a Gitlab account, https://eicweb.phy.anl.gov/users/sign_up, and someone from ATHENA simulation team will approve your account. This account will enable you to copy (git clone) the repo. 

If you do have all items above and still want to try, please continue.

Go to your working dir, e.g., /gpfs02/eic/YOUR_NAME/ATHENA/

## Getting the container

`curl -L get.athena-eic.org | bash`
`./eic-shell`

this command put you in the "container".

## Getting the repo

`git clone git@eicweb.phy.anl.gov:EIC/benchmarks/physics_benchmarks.git`

- switch to branch "diffractive-phi-benchmarks"
`git checkout diffractive-phi-benchmarks`

- local setup
`export JUGGLER_INSTALL_PREFIX=$HOME/stow/juggler # if developing algorithms
export JUGGLER_DETECTOR=athena   # athena is the default
export BEAMLINE_CONFIG=ip6       # ip6 is the default
`
`cd physics_benchmarks
git clone https://eicweb.phy.anl.gov/EIC/benchmarks/common_bench.git setup
source setup/bin/env.sh && ./setup/bin/install_common.sh
source .local/bin/env.sh && build_detector.sh
mkdir_local_data_link sim_output
mkdir -p results
mkdir -p config
`

- temporary needs a calibration:
`bash bin/get_calibrations`

## GEN-SIM-DIGI/RECO-Analysis



This branch you are reading - "diffractive-phi-benchmarks" under physics_benchmark, is currently being developed by Kong.Tu (kongtu@bnl.gov). 

Here the instructions are ONLY intended for people who wants to give it a shot at full ATHENA detector simulation, especially if you do not have much experience in ATHENA (yet). In addition, if you want to follow this particular setup/framework, you need to have ALL of the following:

- BNL RACF account, if you do not have and would like to have one, follow this link, https://wiki.bnl.gov/eic/index.php/Getting_Started. When you just get a new account, you need to setup your EIC environment. If you are desperate enough to use my setup (my setup changes all the time due to other tasks), see /gpfs02/eic/ztu/.cshrc. But I strongly suggest you contact Kolja Kauder <kkauder@gmail.com> (don't tell him you got it from here:) ) 


- Your own ready-to-use HEPMC MC sample.

- Register for a Gitlab account, https://eicweb.phy.anl.gov/users/sign_up, and someone from ATHENA simulation team will approve your account. This account will enable you to copy (git clone) the repo. 

If you do have all items above and still want to try, please continue.
Go to your working dir, e.g., /gpfs02/eic/YOUR_NAME/ATHENA/

## Getting the container

`curl -L get.athena-eic.org | bash`

`./eic-shell`

this command put you in the "container".

## Getting the repo and setups

`git clone git@eicweb.phy.anl.gov:EIC/benchmarks/physics_benchmarks.git`

switch to branch "diffractive-phi-benchmarks"; however, you should create your own branch, own folder, so later can be merged to the master. Mine was created by Wouter Deconinck. Ask him if you are not sure.

`git checkout diffractive-phi-benchmarks`


Local setup (similar to setup codes two levels up.) To be safe, I redo the following setup every time I login, which for some are not necessary obviously.

```
export JUGGLER_INSTALL_PREFIX=$HOME/stow/juggler # if developing algorithms
export JUGGLER_DETECTOR=athena   # athena is the default
export BEAMLINE_CONFIG=ip6       # ip6 is the default
export JUGGLER_N_EVENTS=5 		# number of events to run
```

more setup,

```
cd physics_benchmarks
git clone https://eicweb.phy.anl.gov/EIC/benchmarks/common_bench.git setup
source setup/bin/env.sh && ./setup/bin/install_common.sh
source .local/bin/env.sh && build_detector.sh
mkdir_local_data_link sim_output
mkdir -p results
mkdir -p config
```

(temporary) need a calibration for reco:

`bash bin/get_calibrations`

## Performing GEN-SIM-DIGI/RECO-Analysis steps

- Gen step:

In this example, I have my ready-to-use MC hepmc file, so no need to generate MC. I made this step simply for just doing a copy. Do the following:

`bash benchmarks/diffractive_phi/gen.sh --ebeam 18 --pbeam 110 --config barrel `

You can see from gen.sh, line 73:

`cp /gpfs02/eic/ztu/ATHENA/detectorSimulations/BeAGLE/hepmc3_test_ep_Oct_14/ep_vm.hepmc ${TMP_PATH}/${GEN_TAG}.hepmc`

Besides those printout, you only need to replace this line to copy the hepmc file to `${TMP_PATH}/${GEN_TAG}.hepmc`

- SIM-DIGI-RECO step:

If the first step was done correctly, then this step should work out of the box. Currently ep works, while eA still has issue. 

`bash benchmarks/diffractive_phi/reco_local.sh --ebeam 18 --pbeam 110 --config barrel `

- Analysis step:

This analysis code lives in `benchmarks/diffractive_phi/analysis/diffractive_phi_analysis.cxx` and to run it, simply do:

`bash benchmarks/diffractive_phi/analysis-only.sh --ebeam 18 --pbeam 110 --config barrel`

The output root file will be under `results/diffractive_phi/18on110/barrel_output.root`
Note "barrel" is the config one can choose since the beginning.
