![TREKIS_logo_small](https://github.com/N-Medvedev/TREKIS-3/assets/104917286/e3c4a63a-5b85-497b-93f2-f9a8f17e9bb7)

 ## TREKIS : <ins>T</ins>ime <ins>R</ins>esolved <ins>E</ins>lectron <ins>K</ins>inetics in SHI <ins>I</ins>rradiated <ins>S</ins>olids
 ### Monte-Carlo code modelling electronic kinetics after swift-heavy ion impact on matter
 ### Current version: 3.0.8 (update 17.06.2023)
 
 For all details and instruction, address the files
 !READ_ME_TREKIS_3.doc  or  !READ_ME_TREKIS_3.pdf
 
 ## Disclaimer

Although we endeavour to ensure that the code `TREKIS-3` and results delivered are correct, no warranty is given as to its accuracy. We assume no responsibility for possible errors or omissions. We shall not be liable for any damage arising from the use of this code or its parts or any results produced with it, or from any action or decision taken as a result of using this code or any related material.

This code is distributed as is for non-commercial peaceful purposes only, such as research and education. It is explicitly prohibited to use the code, its parts, its results or any related material for military-related and other than peaceful purposes.
By using this code or its materials, you agree with these terms and conditions.

 ## How to cite

The use of the code is at your own risk. Should you choose to use it, appropriate citations are mandatory:
> 1)	N. A. Medvedev, R. A. Rymzhanov, A. E. Volkov, J. Phys. D. Appl. Phys. 48 (2015) 355303
> 2)	R. A. Rymzhanov, N. A. Medvedev, A. E. Volkov, Nucl. Instrum. Methods B 388 (2016) 41

Should you use this code to create initial conditions for further molecular dynamics simulations of atomic response to the electronic excitation by a swift heavy ion (e.g. with LAMMPS), the following citation is required:

> 3)	R. Rymzhanov, N. A. Medvedev, A. E. Volkov, J. Phys. D. Appl. Phys. 50 (2017) 475301

In a publication, we recommend that at least the following parameters should be mentioned for reproducibility of the results: material, its structure, density, speed of sound, the used CDF coefficients, which processes were included (active) in the simulation, ion type, its energy, the model for SHI charge, number of MC iterations.

## Relevant references
 
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


### 4) Current state of the project, TREKIS-4:

a)	Demonstration of the importance of the nonthermal effects in SHI tracks [11]

b)	Detailed description of the combination of MC-MD, proposed simple model for description of nonthermal effects [12]


## References 

> [1]	N.A. Medvedev, A.E. Volkov, N.S. Shcheblanov, B. Rethfeld, Early stage of the electron kinetics in swift heavy ion tracks in dielectrics, Phys. Rev. B. 82 (2010) 125425. https://doi.org/10.1103/PhysRevB.82.125425
>
> [2]	N.A. Medvedev, A.E. Volkov, B. Rethfeld, N.S. Shcheblanov, Effect of inter-atomic Auger processes on relaxation of electronic vacancies at deep levels of highly ionized atoms in swift heavy ion tracks, Nucl. Instrum. Methods B 268 (2010) 2870–2873. https://doi.org/10.1016/j.nimb.2010.03.021
> 
> [3]	N. Medvedev, Modeling ultrafast electronic processes in solids excited by femtosecond VUV-XUV laser Pulse, AIP Conf. Proc. 582 (2012) 582–592. https://doi.org/10.1063/1.4739911
> 
> [4]	N.A. Medvedev, R.A. Rymzhanov, A.E. Volkov, Time-resolved electron kinetics in swift heavy ion irradiated solids, J. Phys. D. Appl. Phys. 48 (2015) 355303. https://doi.org/10.1088/0022-3727/48/35/355303
> 
> [5]	R.A. Rymzhanov, N.A. Medvedev, A.E. Volkov, Effect of valence holes kinetics on material excitation in tracks of swift heavy ions, Nucl. Instrum. Methods B 365 (2015) 462–467. https://doi.org/10.1016/j.nimb.2015.08.043
> 
> [6]	R.A. Rymzhanov, N.A. Medvedev, A.E. Volkov, Effects of model approximations for electron, hole, and photon transport in swift heavy ion tracks, Nucl. Instrum. Methods B 388 (2016) 41–52. https://doi.org/10.1016/j.nimb.2016.11.002
> 
> [7]	R. Rymzhanov, N.A. Medvedev, A.E. Volkov, Damage threshold and structure of swift heavy ion tracks in Al2O3, J. Phys. D. Appl. Phys. 50 (2017) 475301. https://doi.org/10.1088/1361-6463/aa8ff5
> 
> [8]	R.A. Rymzhanov, N. Medvedev, A.E. Volkov, J.H. O’Connell, V.A. Skuratov, Overlap of swift heavy ion tracks in Al2O3, Nucl. Instrum. Methods B 435 (2018) 121–125. https://doi.org/10.1016/j.nimb.2017.11.014
> 
> [9]	R.A. Rymzhanov, S.A. Gorbunov, N. Medvedev, A.E. Volkov, Damage along swift heavy ion trajectory, Nucl. Instrum. Methods B 440 (2019). https://doi.org/10.1016/j.nimb.2018.11.034
> 
> [10]	R.A. Rymzhanov, N. Medvedev, J.H. O’Connell, V.A. Skuratov, A. Janse van Vuuren, S.A. Gorbunov, A.E. Volkov, Insights into different stages of formation of swift heavy ion tracks, Nucl. Instrum. Methods B 473 (2020) 27–42. https://doi.org/10.1016/j.nimb.2020.04.005
> 
> [11]	N. Medvedev, A.E. Volkov, Nonthermal acceleration of atoms as a mechanism of fast lattice heating in ion tracks editors-pick
Journal of Applied Physics 131 (2022) 225903. https://doi.org/10.1063/5.0095724
> 
> [12]	N. Medvedev, F. Akhmetov, R.A. Rymzhanov, R. Voronkov, A.E. Volkov, Modeling time-resolved kinetics in solids induced by extreme electronic excitation, Adv. Theory Simul. 5, 2200091 (2022). https://doi.org/10.1002/adts.202200091
