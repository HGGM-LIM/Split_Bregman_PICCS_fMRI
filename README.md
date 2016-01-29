# Split_Bregman_PICCS_fMRI
Split Bregman Prior Image-Based Constrained Compressed Sensing (PICCS) for fMRI preclinical data

This repository contains a demo that shows how to use PICCS, which is efficiently implemented with the Split Bregman formulation, for preclinical fMRI, as used in the paper: 

**Chavarrias, C., Abascal, J.F., Montesinos, P. and Desco, M.,"Exploitation of temporal redundancy in compressed sensing reconstruction of fMRI studies with a prior-based algorithm (PICCS)". Med Phys,  42(7): p. 3814 (2015).** 
DOI:% http://dx.doi.org/10.1118/1.4921365

The Split Bregman method separates L2- and L1-norm functionals in such a way that they can be solved analytically in two alternating steps. In the first step a linear system is efficiently solved in the Fourier domain, which can be done in MRI and image denoising problems where operators have representation in the Fourier domain. The computational cost is three FFT per iteration. 

The demo uses cardiac cine small-animal data to simulate an undersampling pattern based on a variable density pdf and compare Spatial TV with Spatiotemporal TV. Both methods are efficiently solved with a computational cost of three FFT per iteration. 

![Statistical map](https://github.com/HGGM-LIM/Split_Bregman_PICCS_fMRI/blob/master/map_full_p_0_01_k_12.jpg)
![Time course at max](https://github.com/HGGM-LIM/Split_Bregman_PICCS_fMRI/blob/master/timecourse.jpg)

The repository contains the following files:

- **fMRI_64x64x115_image.mat:** Absolute image of the central slice of a rat fMRI dataset (115 time frames, healthy rat)
(Acquired data can be found at http://biig.uc3m.es/fmri_rat_data/)

- **map_full_p_0_01_k_12.tif:** Resulting statistical map from the fully sampled dataset

- **timecourse.tif:** The time course of the pixel with highest statistical significance

- **Demo_SpatioTemporalTV_SplitBregman_Sim.m:** This demo loads image, simulate an undersampling pattern for a given acceleration and reconstruct data using spatial total variation (S-TV) and spatiotemporal TV (ST-TV)

- **PICCS_fMRI.m:** PICCS solved using the Split Bregman formulation

- **genRAll.m:** Loop to generate a different sampling patterns for each time point, so that sampling is also pseudo-random along time

- **genPDF.m:** Function to generate a pdf with polynomial variable density sampling, by Michael Lustig 2007, downloaded from (SparseMRI V0.2), http://web.stanford.edu/~mlustig/SparseMRI.html, see M. Lustig, D.L

- **genSampling_LIM.m:** Monte-carlo algorithm to generate a sampling pattern. Modified from the original function by Michael Lustig 2007

- **maxSidePeakRatio.m:** Computes the maximum sidelobe to peak ratio of point spread function for undersampled data. Used within genSampling_LIM.m


If you need to contact the author, please do so at cristina.chavarrias@gmail.com, juanabascal78@gmail.com, desco@hggm.es
