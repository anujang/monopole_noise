# monopole_noise
Modelling the static and dynamic change of monopole movement in a 2D square lattice, with changing temperature and applied magnetic field,

In this repo, there are jupyter notebook files, which can be run to find a variety of things, such as the behaviour of monopoles in 2D lattices. Due to noise variation, ideally, this is run on a HCP and averaged using the bash script provided. The averaged data file is then exported from the HCP and analysed using the python files to obtain useful graphs, e.g., change in power spectral density aka noise. 

A monte carlo method is used to simulate the behaviour of the monopoles utilising the Hamilonian of systems including the Ising and nearest neighbour models. 

This result was later then independantly confirmed by this paper https://arxiv.org/pdf/2211.09784.pdf. //
This code and results were part of my 4th year dissertation for my Physics MPhys at Cardiff University. Initially I had planned to publish these results but this paper had beat us to it. 
