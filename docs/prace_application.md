# Some thoughts about a future PRACE time allocation

The application for support with code optimization should be written with a future PRACE time allocation in mind. What should the code be able to do? And how much could be accomplished within a typical PRACE projec? 

## Hard numbers:

* All numbers refer to model sizes on the order of 512^3 
* Time used with current RT code version for 1 wavelength on 1 core: 3 hours 
* Size of typical PRACE allocation: 10 million (core) hours 
* Anticipated RT code speed-up: 4 

## Soft numbers: 

* ALMA - Wavelengths in one receiver band: 4000 
* ALMA - bands to be considered: 2 - 4 (out of eventually 10) 
* Bifrost - duration of simulation: 1800 s (- 3600 s) 
* Bifrost - output cadence (= RT input): 1 s (-2 s)
* Bifrost - number of different models to be considered: 6 - 10

## Resulting computation time 

* 1 Bifrost model, 1800s with 1s cadence, 2 bands with each 4000 channels: 1 x 1800 x 4000 = 14.4 million
* Total time on a single core (in million hours) with current RT code = 43.2 million hours 
* Total time on a single core (in million hours) with optimized RT code = 10.8 million hours 

## Conclusion

With the current numbers, a typical PRACE allocation will only allow for 2 ALMA bands for 1 Bifrost model. 

===

# Context
	Prace application for _Code development, starting on Tier-1 system with support from PRACE experts_

# General informations
	**Start date:  01.05.2017** (_The start date should be either one and half or two months after the announced date for evaluation cut-off._)

# Requested computer systems
	It will be *Tier-1 Abel* for us. [Abel](https://www.uio.no/english/services/it/research/hpc/abel/more/)

## 1. Summary of the project 
(To be published in the PRACE website. Maximum 500 words)

**Version below has 488 words.**

The Atacama Large Millimeter/submillimeter Array (ALMA), which is currently the world's largest astronomical observatory, opened up a new window on the universe. The interferometric array is located on the Chajnantor plateau in the Chilean Andes at an altitude of 5000 m and consists of 66 antennas, most of them with a diameter of 12 m. By combining the antennae, they act like a giant telescope with baselines of up to 16 km. Already within the first few years of operation, ALMA led to many scientific discoveries. Since December 2016, ALMA is also used for observations of our Sun. It observes the Sun at a spatial resolution, which is unprecedented in this wavelength range, and offers novel means of determining the properties of the plasma in the Sun's outer atmospheric layers. Due to the properties of the solar radiation at millimeter wavelengths, ALMA serves as a linear thermometer, mapping narrow layers at different heights. It can measure the thermal structure and dynamics of the solar atmosphere and thus sources and sinks of atmospheric heating.

Among other expected scientific results, ALMA promises significant steps towards understanding the intricate dynamics and physical processes that, in combination, might yield the solution of the coronal heating problem - a long standing fundamental question in modern astrophysics. However, ALMA's novel diagnostic capabilities, which will ultimately progress our understanding of the Sun, still need to be developed and understood further in order to fully exploit their potential. Detailed numerical simulations of the solar atmosphere and artificial observations of the Sun play a key role in this respect. 

Such artificial observations of the Sun will be produced as part of the SolarALMA project at the University of Oslo, which is funded with a Consolidator Grant by the European Research Council (ERC), in collaboration with Dr. de la Cruz Rodriguez from the University of Stockholm. The overall aim of the SolarALMA project is to utilize the first observations of the Sun with ALMA and to develop the required observing and analysis techniques. An important step in this endeavour is the development of realistic numerical models of the solar atmosphere, which can be used to test how to optimally set up solar observations and how to analyse and interpret such observations. While 3D numerical models are routinely produced on high-performance computers, the codes available for producing the corresponding (artificial) ALMA observations of such models did not perform well enough so far. A new code has therefore been developed, which solves the radiative transfer equation for a 3D numerical model and thus reveals how the modeled part of the Sun would look like through the eyes of ALMA at (sub)millimeter wavelengths. The new code is in an advanced stage but still needs to be optimized in order to be capable of carrying out for the foreseen computations and thus providing the basis for essential studies, which will result in optimizing ALMA's scientific impact for observing and understanding the Sun. 

## 2. Scientific case of the project
Explain the scientific case for which you intend to use the code(s). Maximum 500 words.

**Version below has 500 words.**

In order to understand the diagnostic potential of the Atacama Large Millimeter/submillimeter Array (ALMA) for observations of the Sun, realistic models and corresponding artificial observations are needed. In summary, we would use a 3D model of the solar atmosphere as input for a radiative transfer (RT) code, and then apply ALMA's instrumental properties to the RT results. Changing the parameters of the last step (e.g., the temporal resolution of the observation) then allows for investigating the impact of the change on the artificial observations and the extent to which initially contained information can be recovered from the artificial observations. The fundamental strength of our approach is that we can compare the processed artificial observations with the original solar input model. Such comparisons directly reveal, which aspects of the input model can be recovered or not and what is needed to get the most out of the data. The resulting strategies can be applied directly to real ALMA observations of the Sun and yield the potential of greatly boosting ALMA's diagnostic potential, promising high-impact results. In particular, we plan to look for observational signatures of small-scale, intermittent events in the solar atmosphere such as, e.g., various oscillations and wave modes and nano-flares, which could contribute to the still not fully explained heating of the outer layers of our Sun. Finding such small-scale events demands maximizing the capabilities of ALMA's solar observing modes. This task can be substantially supported by the numerical simulations described here but requires detailed RT calculations, which currently do not exist.  

Suitable 3D models of the solar atmosphere have been produced with the Bifrost code in Oslo and will serve as input for a new radiative transfer code that calculates how the model looks like at the wavelengths observed with ALMA. As a first step, we plan to simulate an ALMA observation run in all detail. That includes observations in two ALMA receiver bands (i.e., frequency/wavelength ranges), which each can have around 4000 spectral channels. ALMA observes at very high cadence, which will be 1s in the cycle starting in October 2017. We plan therefore to process a time series from the 3D Bifrost model with 1s cadence over a period of 30min. In total, the radiative transfer code has to compute the emergent intensity at 2x4000 wavelengths for 1800 time steps. This task is computationally expensive and practically impossible with previously available radiative transfer codes that required substantial human interaction. The new code developed for this task has the potential to handle the challenging computations, given that it can be sped up by a factor of ideally 4, which we consider perfectly possible.  

Among other things, we plan to address the following points with the anticipated data set: 1) Tests and improvements of ALMA data calibration procedures. 2) Multi-frequency synthesis (boosting spatial resolution by sacrificing spectral resolution). 3) Novel "burst mode" (combining consecutive time steps for higher image fidelity). 

Finally, it should be noted that the same tools can be applied to models of other stars, too. 

## 3. Computer resources requested  (for Types A, B, C), or expected for Type D
The input for our radiative transfer (RT) code are 3D solar/stellar atmosphere models (mainly those generated with the Bifrost code), which are sets of 3D arrays of floats for typically 6 primary physical variables (most of the time in single precision). For a typical Bifrost model, the size of the input accounts to 768^3 * 6 * 64bits = ~ 22GB for one snapshot. The RT output, in extreme cases will be of the same order. Thus in total we need at least 1TB of storage for input & output, to be able test our code with short time series with a few snapshots.  
	
Total storage required (Gbyte) : 1014 (1TB)  
Maximum amount of memory per core (Mbyte)	**~0.5GB**

## 4. Please provide the details listed below for the main simulation application  (for Types A, B, C, D)

	Name and version:  **RTam** 
	         Webpage:  **none**

It is currently an alpha version without license assigned but we plan to release it under GPL or MIT license.

## 5. Describe the main algorithms and how they have been implemented and parallelized  

At the core of computation is a nonlinear solver for radiative transfer (RT), which requires the solution of the (polarized) radiative transfer problem. This problem involves the solution of a set of coupled first-order linear differential equations, one for each Stokes parameter.

Given that each column in the atmosphere model is a completely independent calculation, the problem is very suitable for parallel computing. A very limited amount of communication is necessary between different processes. The code has been parallelized using MPI. The total number of columns to be calculated (tasks) is divided by the number of processes, so that the amount of tasks each process has is about the same and is known in the beginning of the execution. Each process starts working through its task list until it is finished, and then writes the output to the disk and waits for the other processes to complete.

## 6. Current and target performance  (for Types A, B, C, D ; including the points below. Maximum 250 words)
The code is naively parallel but it is very computationally demanding. It scales with ease up to 512 cores but it has not been tested yet  with a higher number of CPUs. We believe that better performance can be achieved with a hybrid approach (OpenMP/MPI) or by using accelerators like Intel Xeon Phi or GPGPU. 
	
What is the target for scalability and performance? (i.e. what performance is needed to reach the envisaged scientific goals) 

While the scalability shouldn't be a problem we would need a speed up by a factor of ideally 4 to be able generate the artificial observations for ALMA as described above (item 2). 

## 7. Confidentiality  (for Types A, B, C, D)
Is any part of the project covered by confidentiality?    No

## 8. Describe the I/O strategy regarding the parameters indicated below  (for Types A, B, C, D)
a) Is I/O expected to be a bottleneck? (Maximum 50 words) 
Although the input and output is very volumetric, the time spend on I/O operations is negligible when compared to the computation time. 
b) Implementation: I/O libraries, MPI I/O, netcdf, HDF5 or other approaches. (Maximum 50 words) We use the parallel HDF5 library for I/O.
c) Frequency and size of data output and input (Maximum 50 words) Input read only at the beginning (~25GB for snapshot),  output saved only at the end of computation (~20GB for each snapshot)
d) Number of files and size of each file in a typical production run (Maximum 50 words) 1 hdf5 file (~25GB -> xTB) for input and 1 for output (~20GB -> xTB). The code can process time series a model snapshots sequentially. The large production run mentioned under item 2 would thus result in 1800 I/O files (one at a time).  

## 9. Main performance bottlenecks  (for Types B, C, D. Maximum 250 words)
Currently the slowest part is the EOS (Equation of State) module which strongly depends on expensive calculation of integrals of exponential functions.

## 10. Describe possible solutions you have considered to improve the performance of the project   
	(for Types B, C, D. Maximum 250 words)
The major improvement would require changes in the algorithm. We believe that using interpolation tables for the equation of state instead of the explicit calculations of integrals of exponential functions would speed up the code substantially keeping reasonable precision of the solution.  

The code hasn't been optimized so far for better vectorization, thus in depth performance analysis is needed with fine tuning (loop by loop) for better SIMD usage. We also believe that our code should be easy to accelerate on Xeon Phi (KNL) or/and GPGPU.  

Adding OpenMP should also reduce memory overhead and would improve scalability. 

## 11. Describe the application enabling/optimization work that needs to be performed to achieve the target performance  (for Types B, C, D. Maximum 250 words)
In depth performance analysis. 
Vectorization guidance (e.g., via OpenMP 4.0) 
Code adaptation to GPGPU/KNL 
Introducing OpenMP

## 12. Which computational performance limitations do you wish to solve with this project?   
	(for Types B, C, D. Maximum 250 words)
Reduce number of memory bound computations.
Improve vectorization (increase FLOPs)
Validate whether it would be possible to port this code to GPGPU 
Performance estimate for/on Intel Xeon Phi KNL 

## 13. Describe the impact of the optimization work proposed?  (for Types B, C, D. Maximum 250 words)
	* Is the code widely used?

	* Would the code be used only within this original research project?
No. Many applications are possible and foreseen, not only for modelling of the Sun but also other stars. 

* Would the code be used for other similar research projects with minor modifications?
	
Yes, as part of the ongoing SolarALMA project, stellar simulations and also in collaboration with external international collaborators.  	
	
	* Would the code be used in many research projects of the research field indicated in the proposal?
	
	Yes. We intend to make the code available for a large user base, in particular for the Solar Simulations for the Atacama Large Millimeter Observatory Network (SSALMON), which hosted and organized by our group. It currently involves 86 scientists from around the world.  
	
	* Would the modification be easy to add to the main release of the software?
	Yes. 
	

## 14. Describe the request plans for work with support from PRACE experts  (for Type C, D)
	

### a. Describe the level of collaboration with PRACE experts you have planned for and how much effort (person months) have you reserved for this?

	Specify a rough estimate for the amount of person months this work entails :  4PMs

### b. Describe the optimization work you expect to be done with the support of a PRACE expert for your project
	While we have resources to work on improving the algorithm and introducing OpenMP, we would like a PRACE expert to help us with performance analysis and vectorization. We also would like to work in sync on porting our code to GPGPU/KNL. 

### c. Please specify the amount of PRACE experts
	person months required to support your project (1-6 PMs):  4MPs

## 15. Please describe your hardware architecture currently in use to allow an optimal selection of a feasible Tier-1 system  (for Type D only, Maximum 250 words)
	 We regularly run our codes on multi-core, Intel Xeon, 64GB/compute node
	 Some of them are equipped with Infiniband. 
	 We are planning on buying a small cluster with GPU (2x Nvidia p100) and Xeon KNL. Access to a similar system where we could test the adapted/optimized code beforehand would therefore be appreciated. 

 
