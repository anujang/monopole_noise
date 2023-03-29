# Monopole Noise In Artificial Spin Ices
Modelling the static and dynamic change of monopole movement in a 2D square lattice, with changing temperature and applied magnetic field,

In this repo, there are jupyter notebook files, which can be run to find a variety of things, such as the behaviour of monopoles in 2D lattices. Due to noise variation, ideally, this is run on a HCP and averaged using the bash script provided. The averaged data file is then exported from the HCP and analysed using the python files to obtain useful graphs, e.g., change in power spectral density aka noise. 

A monte carlo method is used to simulate the behaviour of the monopoles utilising the Hamilonian of systems including the Ising and nearest neighbour models. The material simulated is called a spin ice, which exists in real life, though very expensive. So it is simulated in 2D, since the effects leading to 'magnetic monopoles' can be seen clearly. 

The aim of this project was to characterise the behaviour of noise spectral density of 2D nearest-neighbour spin ice model. This was to identify which of Kirschner et al [1] or Samarakoon et al [2] were correct in their claims for the noise spectral density. 

# Nearest Neighbour Modal
This is the main model used, where atomic spins live along the edge of the vertices, as shown in the picture. 

<img width="287" alt="Screenshot 2023-03-29 at 18 32 33" src="https://user-images.githubusercontent.com/99356066/228620808-b017b0d0-5aac-4551-9d83-7f77cc9102c2.png">

This modal can be easily used to there being a finite number of configurations possible for the monopole orientations. This is known as the 16 vertex model, in the following diagram. 

<img width="433" alt="Screenshot 2023-03-29 at 18 33 41" src="https://user-images.githubusercontent.com/99356066/228621086-3ccde1f9-ad2f-41d7-aaaa-4f4d2e6151d8.png">
[D. Levis. Two-dimensional Spin Ice and the Sixteen-Vertex Model. PhD thesis, Universite Pierre et Marie Curie - Paris VI, (2012).]

# Monte Carlo Method
Since all configurations are known, a single point is initally chosen and its energy is found. Then a flip of spin is attempted and the cost of energy is found, only if the cost is negative, the flip is allowed. (for more info, see disseration/powerpoint in the dissertation folder)
Each monte carlo step reppresents one time step in real life. To compare the two papers mentioned at the start, the measurement of total magnetisation in the lattice is measured once each MC time step and each time a flip is attempted. Bearing in mind, in each MC step, L^2 flips are attempted (where L is length of lattice)

<img width="281" alt="Screenshot 2023-03-29 at 18 40 14" src="https://user-images.githubusercontent.com/99356066/228622545-53dfd5eb-e421-4144-b394-b1f2c1162d6c.png">

[M. Goryca, X. Zhang, J. Li, A. L. Balk, J. D. Watts, C. Leighton, C. Nisoli, P. Schiffer, and S. A. Crooker. Field-induced magnetic monopole plasma in artificial spin ice. Phys. Rev. X, 11:011042, March 2021.]

The dynamics of the monopoles are found by apply a magnetic field onto the material and varying its strenght, which is emulated via $|B_x| + |B_y| \propto 4(J_1 - J_2) = 3.2K_b T $

# Power Spectral Density
Noise in system found through Fourier transform of autocorrelation function:

𝐶(𝜏)= ⟨𝑀(𝑡)𝑀(𝑡+𝜏)⟩

𝑃(𝜔)=𝐶 ̂(𝜔)=1/𝐿 |∫_0^𝐿〖𝑀(𝑡) 𝑒^𝑖𝜔𝑡 𝑑𝑡〗|^2

From this the decay exponent is found to be simply

P(\omega) \propto 1/\omega^\beta

Common decay exponent/noise values are known to be:

0 - White

1 - Pink

2 - Red/Brownian, i.e., random motion

## Nyquist frquency

This is an important topic to understand for this project

Arises from the rate at which we must sample for Fourier transform 

Aliasing: undersampling of measurements leading to incorrect reconstruction of data.

Too slow, get aliasing as seen here

<img width="405" alt="Screenshot 2023-03-29 at 18 46 43" src="https://user-images.githubusercontent.com/99356066/228624212-d0501c01-28ce-4fc3-90e0-808164df948c.png">


Too fast, oversample and go past real Nyquist frequency and get nonsense results.
Past the Nyquist frequency, data reflects due to Fourier transform..

<img width="404" alt="Screenshot 2023-03-29 at 18 47 13" src="https://user-images.githubusercontent.com/99356066/228624329-94edafd4-9e29-49d5-bf1d-7422aef2f8b0.png">

# Measurement per MC Step

Measuring once per MC step gave a noise of $ \beta = 0.950 +- 0.017 $, which agrees with literature, and indicates pink noise, which makes sense for the probabilistic nature of magnetic monopoles.

<img width="508" alt="image" src="https://user-images.githubusercontent.com/99356066/228624480-f895f2bc-3888-4509-9896-3428f5aec5bb.png">


# Measurement per Attempted Flip

Before Nyquist Frequency: 𝛽=0.947±0.019
After Nyquist Frequency: 𝛽=1.965±0.021

As seen from the results here, when one measures once per attempted flip, there is an overmeasurement by a factor of $L^2$. This is seen as if you take the plateau in the same place as before and measure the beta value, it is found to be around 2, which makes no sense for a system of magentic monopoles. If the Nyquist frequecy is remembered and applied, it is obvious there is oversampling by this factor, and finding the decay exponent before the nonsensical data past the Nyquist frequency, yields values of around 1, which seems to agree with the results obtained in the previous section. This is further seen in the follow section, across a range of applied magnetic fields (as a function of temperature)


<img width="526" alt="image" src="https://user-images.githubusercontent.com/99356066/228624623-52b97f19-6f3b-496e-9976-1bc801cd6480.png">

#Comparison of Methods

It becomes more obvious in these images, that the method from Samarakoon et al, was infact wrong, they were only measuring the mathematical nonsense which occurs past the Nyquist frequency, which simply occurs as a byproduct of the fourier transform invovled in the MC steps. The true decay constant nature of a monopole system lies around 1. It is well established in literature that monopole systems are not random (i.e., an unconstrained random walk system will have a decay exponent of 2, as it is entirely random). Monopoles have memory and therefore this decay exponent must always be below 2. 

<img width="644" alt="image" src="https://user-images.githubusercontent.com/99356066/228624789-0115939e-ed17-470f-96d0-94c57b84f5bb.png">


<img width="640" alt="image" src="https://user-images.githubusercontent.com/99356066/228624821-1050b6a1-ce4e-4b4b-89b7-42fd1a26ae35.png">


This result was later then independantly confirmed by this paper https://arxiv.org/pdf/2211.09784.pdf. 

This code and results were part of my 4th year dissertation for my Physics MPhys at Cardiff University. Initially I had planned to publish these results but this paper had beat us to it. 




[1] F. K. K. Kirschner et al. Proposal for the detection of magnetic monopoles in spin ice via nanoscale magnetometry. Phys. Rev. B, 97:140402, (2018).
[2] A. M. Samarakoon et al. Anomalous magnetic noise in an imperfectly flat landscape in the topological magnet Dy2Ti2O7. Proceedings of the National Academy of Sciences, 119(5), (2022).
