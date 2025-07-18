![TREKIS_logo_small](https://github.com/N-Medvedev/TREKIS-3/assets/104917286/e3c4a63a-5b85-497b-93f2-f9a8f17e9bb7)

 ## TREKIS : <ins>T</ins>ime <ins>R</ins>esolved <ins>E</ins>lectron <ins>K</ins>inetics in SHI <ins>I</ins>rradiated <ins>S</ins>olids
 [![DOI](https://zenodo.org/badge/490195185.svg)](https://zenodo.org/badge/latestdoi/490195185)
 
 ### Monte-Carlo code modelling electronic kinetics after swift-heavy ion impact on matter
 ### Current version: 3.3.0 (update 18.07.2025)
 
 For all details and instruction, address the files
 !READ_ME_TREKIS_3.doc  or  !READ_ME_TREKIS_3.pdf
 
 ## Disclaimer

_This code is work in progress, anything might change without a notice, bugfixes and patches are expected!_

Although we endeavour to ensure that the code TREKIS-3 and results delivered are correct, no warranty is given as to its accuracy (for details, see GPL-3.0 license). This code was developed for non-commercial peaceful purposes only, such as research and education.

 ## How to cite

The use of the code is at your own risk. Should you choose to use it, please cite the appropriate works:
> 1)	N. A. Medvedev, R. A. Rymzhanov, A. E. Volkov, J. Phys. D. Appl. Phys. 48 (2015) 355303
> 2)	R. A. Rymzhanov, N. A. Medvedev, A. E. Volkov, Nucl. Instrum. Methods B 388 (2016) 41

Should you use this code to create initial conditions for further molecular dynamics simulations of atomic response to the electronic excitation by a swift heavy ion (e.g. with LAMMPS), please also cite the following work:

> 3)	R. Rymzhanov, N. A. Medvedev, A. E. Volkov, J. Phys. D. Appl. Phys. 50 (2017) 475301

In a publication, we recommend that at least the following parameters should be mentioned for reproducibility of the results: material, its structure, density, speed of sound, the used CDF coefficients, which processes were included (active) in the simulation, ion type, its energy, the model for SHI charge, number of MC iterations.

## Relevant references

### 0) Review of modeling of swift-heavy ion effects in matter:

N. Medvedev, A.E. Volkov, R. Rymzhanov, F. Akhmetov, S. Gorbunov, R. Voronkov, P. Babaev, _Frontiers, challenges, and solutions in modeling of swift heavy ion effects in materials_,
J. Appl. Phys. 133, 100701 (2023)
https://doi.org/10.1063/5.0128774

 
### 1) Early and preliminary works, prior to TREKIS:

a) Monte Carlo simulations of selected cases, demonstrating importance of nonequilibrium electron kinetics in SHI tracks [1,2]

b) Methodology of CDF cross sections [3]


### 2) Presentation of TREKIS-3:

a) Detailed description of the code [4,5]

b) Analysis of model parameters (and first added photon transport) [6]


### 3) Combination of TREKIS with MD simulations (LAMMPS):

a)	First implementation and validation of the methodology [7,8]

b)	Description of the damage along entire SHI track, no fitting parameters [9]

c)	A brief review of the TREKIS+MD results [10]


### 4) [TREKIS-4](https://github.com/N-Medvedev/TREKIS-4):

_Note that TREKIS-4 is a separate project, developed independently and in parallel with TREKIS-3_

a)	Demonstration of the importance of the nonthermal effects in SHI tracks [11]

b)	Detailed description of the combination of MC-MD, proposed simple model for description of nonthermal effects [12]


## References 

> [1]	N.A. Medvedev, A.E. Volkov, N.S. Shcheblanov, B. Rethfeld, _Early stage of the electron kinetics in swift heavy ion tracks in dielectrics_, Phys. Rev. B. 82 (2010) 125425. https://doi.org/10.1103/PhysRevB.82.125425
>
> [2]	N.A. Medvedev, A.E. Volkov, B. Rethfeld, N.S. Shcheblanov, _Effect of inter-atomic Auger processes on relaxation of electronic vacancies at deep levels of highly ionized atoms in swift heavy ion tracks_, Nucl. Instrum. Methods B 268 (2010) 2870–2873. https://doi.org/10.1016/j.nimb.2010.03.021
> 
> [3]	N. Medvedev, _Modeling ultrafast electronic processes in solids excited by femtosecond VUV-XUV laser Pulse_, AIP Conf. Proc. 582 (2012) 582–592. https://doi.org/10.1063/1.4739911
> 
> [4]	N.A. Medvedev, R.A. Rymzhanov, A.E. Volkov, _Time-resolved electron kinetics in swift heavy ion irradiated solids_, J. Phys. D. Appl. Phys. 48 (2015) 355303. https://doi.org/10.1088/0022-3727/48/35/355303
> 
> [5]	R.A. Rymzhanov, N.A. Medvedev, A.E. Volkov, _Effect of valence holes kinetics on material excitation in tracks of swift heavy ions_, Nucl. Instrum. Methods B 365 (2015) 462–467. https://doi.org/10.1016/j.nimb.2015.08.043
> 
> [6]	R.A. Rymzhanov, N.A. Medvedev, A.E. Volkov, _Effects of model approximations for electron, hole, and photon transport in swift heavy ion tracks_, Nucl. Instrum. Methods B 388 (2016) 41–52. https://doi.org/10.1016/j.nimb.2016.11.002
> 
> [7]	R. Rymzhanov, N.A. Medvedev, A.E. Volkov, _Damage threshold and structure of swift heavy ion tracks in Al2O3_, J. Phys. D. Appl. Phys. 50 (2017) 475301. https://doi.org/10.1088/1361-6463/aa8ff5
> 
> [8]	R.A. Rymzhanov, N. Medvedev, A.E. Volkov, J.H. O’Connell, V.A. Skuratov, _Overlap of swift heavy ion tracks in Al2O3_, Nucl. Instrum. Methods B 435 (2018) 121–125. https://doi.org/10.1016/j.nimb.2017.11.014
> 
> [9]	R.A. Rymzhanov, S.A. Gorbunov, N. Medvedev, A.E. Volkov, _Damage along swift heavy ion trajectory_, Nucl. Instrum. Methods B 440 (2019). https://doi.org/10.1016/j.nimb.2018.11.034
> 
> [10]	R.A. Rymzhanov, N. Medvedev, J.H. O’Connell, V.A. Skuratov, A. Janse van Vuuren, S.A. Gorbunov, A.E. Volkov, _Insights into different stages of formation of swift heavy ion tracks_, Nucl. Instrum. Methods B 473 (2020) 27–42. https://doi.org/10.1016/j.nimb.2020.04.005
> 
> [11]	N. Medvedev, A.E. Volkov, _Nonthermal acceleration of atoms as a mechanism of fast lattice heating in ion tracks_, Journal of Applied Physics 131 (2022) 225903. https://doi.org/10.1063/5.0095724
> 
> [12]	N. Medvedev, F. Akhmetov, R.A. Rymzhanov, R. Voronkov, A.E. Volkov, _Modeling time-resolved kinetics in solids induced by extreme electronic excitation_, Adv. Theory Simul. 5, 2200091 (2022). https://doi.org/10.1002/adts.202200091
