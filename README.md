# TREKIS-3
 TREKIS-3: Time-Resolved Electron Kinetics in SHI-Irradiated Solids -- a Monte Carlo simulation of the effects of SHI impact on matter
 
 For all details and instruction, address the files
 !READ_ME_TREKIS_3.doc  or  !READ_ME_TREKIS_3.pdf
 
 *Disclaimer*

Although we endeavour to ensure that the code TREKIS-3 and results delivered are correct, no warranty is given as to its accuracy. We assume no responsibility for possible errors or omissions. We shall not be liable for any damage arising from the use of this code or its parts or any results produced with it, or from any action or decision taken as a result of using this code or any related material.

This code is distributed as is for non-commercial peaceful purposes only, such as research and education. It is explicitly prohibited to use the code, its parts, its results or any related material for military-related and other than peaceful purposes.
By using this code or its materials, you agree with these terms and conditions.

 *How to cite*

The use of the code is at your own risk. Should you choose to use it, appropriate citations are mandatory:
1)	N. A. Medvedev, R. A. Rymzhanov, A. E. Volkov, J. Phys. D. Appl. Phys. 48 (2015) 355303
2)	R. A. Rymzhanov, N. A. Medvedev, A. E. Volkov, Nucl. Instrum. Methods B 388 (2016) 41

Should you use this code to create initial conditions for further molecular dynamics simulations of atomic response to the electronic excitation by a swift heavy ion (e.g. with LAMMPS), the following citation is required:

3)	R. Rymzhanov, N. A. Medvedev, A. E. Volkov, J. Phys. D. Appl. Phys. 50 (2017) 475301

In a publication, we recommend that at least the following parameters should be mentioned for reproducibility of the results: material, its structure, density, speed of sound, the used CDF coefficients, which processes were included (active) in the simulation, ion type, its energy, the model for SHI charge, number of MC iterations.
 
 *Relevant references*
 
*1)* Early and preliminary works, prior to TREKIS:
a) Monte Carlo simulations of selected cases, demonstrating importance of nonequilibrium electron kinetics in SHI tracks [1,2]
b) Methodology of CDF cross sections

*2)* Presentation of TREKIS-3:
a) Detailed description of the code [3,4]
b) Analysis of model parameters (and first added photon transport) [5]

*3)* Combination of TREKIS with MD simulations (LAMMPS):
a) First implementation and validation of the methodology [6,7]
b) Description of the damage along entire SHI track, no fitting parameters [8]
c) A brief review of the TREKIS+MD results [9]

*4)* Current state of the project, TREKIS-4:
a) Demonstration of the importance of the nonthermal effects in SHI tracks [10]
b) Detailed description of the combination of MC-MD, proposed simple model for description of nonthermal effects [11]

 *Bibliography*
 
[1]	N.A. Medvedev, A.E. Volkov, N.S. Shcheblanov, B. Rethfeld, Early stage of the electron kinetics in swift heavy ion tracks in dielectrics, Phys. Rev. B. 82 (2010) 125425. https://doi.org/10.1103/PhysRevB.82.125425

[2]	N.A. Medvedev, A.E. Volkov, B. Rethfeld, N.S. Shcheblanov, Effect of inter-atomic Auger processes on relaxation of electronic vacancies at deep levels of highly ionized atoms in swift heavy ion tracks, Nucl. Instruments Methods Phys. Res. Sect. B Beam Interact. with Mater. Atoms. 268 (2010) 2870–2873. https://doi.org/10.1016/j.nimb.2010.03.021

[3]	N.A. Medvedev, R.A. Rymzhanov, A.E. Volkov, Time-resolved electron kinetics in swift heavy ion irradiated solids, J. Phys. D. Appl. Phys. 48 (2015) 355303. https://doi.org/10.1088/0022-3727/48/35/355303

[4]	R.A. Rymzhanov, N.A. Medvedev, A.E. Volkov, Effect of valence holes kinetics on material excitation in tracks of swift heavy ions, Nucl. Instruments Methods Phys. Res. B. 365 (2015) 462–467. https://doi.org/10.1016/j.nimb.2015.08.043

[5]	R.A. Rymzhanov, N.A. Medvedev, A.E. Volkov, Effects of model approximations for electron, hole, and photon transport in swift heavy ion tracks, Nucl. Instruments Methods Phys. Res. Sect. B Beam Interact. with Mater. Atoms. 388 (2016) 41–52. https://doi.org/10.1016/j.nimb.2016.11.002

[6]	R. Rymzhanov, N.A. Medvedev, A.E. Volkov, Damage threshold and structure of swift heavy ion tracks in Al2O3, J. Phys. D. Appl. Phys. 50 (2017) 475301. https://doi.org/10.1088/1361-6463/aa8ff5

[7]	R.A. Rymzhanov, N. Medvedev, A.E. Volkov, J.H. O’Connell, V.A. Skuratov, Overlap of swift heavy ion tracks in Al2O3, Nucl. Instruments Methods Phys. Res. Sect. B Beam Interact. with Mater. Atoms. 435 (2018) 121–125. https://doi.org/10.1016/j.nimb.2017.11.014

[8]	R.A. Rymzhanov, S.A. Gorbunov, N. Medvedev, A.E. Volkov, Damage along swift heavy ion trajectory, Nucl. Instruments Methods Phys. Res. Sect. B Beam Interact. with Mater. Atoms. 440 (2019). https://doi.org/10.1016/j.nimb.2018.11.034

[9]	R.A. Rymzhanov, N. Medvedev, J.H. O’Connell, V.A. Skuratov, A. Janse van Vuuren, S.A. Gorbunov, A.E. Volkov, Insights into different stages of formation of swift heavy ion tracks, Nucl. Instruments Methods Phys. Res. Sect. B Beam Interact. with Mater. Atoms. 473 (2020) 27–42. https://doi.org/10.1016/j.nimb.2020.04.005

[10]	N. Medvedev, A.E. Volkov, Reconciling anomalously fast heating rate in ion tracks with low electron-phonon coupling, (2021). https://arxiv.org/abs/2109.04401v1

[11]	N. Medvedev, F. Akhmetov, R.A. Rymzhanov, R. Voronkov, A.E. Volkov, Modeling time-resolved kinetics in solids induced by extreme electronic excitation, 2022. https://arxiv.org/abs/2201.08023v1