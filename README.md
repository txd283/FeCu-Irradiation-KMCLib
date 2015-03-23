Welcome to my repository for my third year undergraduate project at The University of Birmingham, UK. This code is for modelling irradiation damage of RPV steels with on-lattice Kinetic Monte Carlo technique.

I have used KMCLib (latest commit 9ce74cf9f27a29bffd8020194e225d36c384a854) https://github.com/leetmaa/KMCLib for the backend KMC.  I have included a complied version for OS X 10.10 64bit (clang no MPI) and Ubuntu 14.04 64bit (gcc no MPI) in KMCLib > Backend > Backend.so.<version>

This is the complete model for my project (23/3/15). It does

- Fe-Cu alloy with vacancies
- Vacancy diffusion with first and second neighbours
- Vacancy clustering
- Copper clustering
- CFG format output support for atomeye visualisation
- Clustering calculator to generate a distribution of cluster size

Planned by never implemented:

- include interstitials and their annihilation with vacancies (nearest neighbours)
- Dislocations (sinks)
- Nickel
- Incorporate molecular dynamics damage cascades

Open-source and under the MIT Licence. If you have any questions, please feel free to email me at txd283@bham.ac.uk
