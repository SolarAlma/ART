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

The Atacama Large Millimeter/submillimeter Array (ALMA), which is currently the worlds' largest astronomical observatory, opened up a new window on the universe. Already within the first few years of operation ALMA led to many scientific discoveries. Since December 2016, ALMA is also used for observations of our Sun. It observes the Sun at a spatial resolution, which is unprecedented in this wavelength range, and offers novel means of determining the properties of the plasma in the Sun's outer atmospheric layers. Among other expected scientific results, ALMA promises significant steps towards understanding the intricate dynamics and physical processes that, in combination, might yield the solution to the coronal heating problem - a long standing fundamental question in modern astrophysics. However, ALMA's novel diagnostic capabilities, which will ultimately progress our understanding of the Sun, still need to be developed and understood further in order to fully exploit their potential. Detailed numerical simulations of the solar atmosphere and artificial observations of the Sun play a key role in this respect. 

Such artificial observations of the Sun will be produced as part of the SolarALMA project, which is funded with a Consolidator Grant by the European Research Council (ERC). The overall aim of the project is to utilize the first observations of the Sun with ALMA and to develop the required observing and analysis techniques. An important step in this endeavour is the development of a realistic numerical models of the solar atmosphere, which can be used to test how to optimally set up solar observations and how to analyse and interpret such observations. While 3D numerical models are routinely produced on high-performance computers, the codes available for the corresponding (artificial) observations of such models did not perform well enough. A new code has therefore been developed, which solves the radiative transfer equation for a 3D numerical model and thus reveals how the modelled part of the Sun would look like through the eyes of ALMA at (sub)millimeter wavelengths. The new code is in advanced stage now but still needs to be optimized in order to be capable of carrying out for the foreseen computations. 


## 2. Scientific case of the project
Explain the scientific case for which you intend to use the code(s). Maximum 500 words.
**Version below at 553 words. Needs to be cut.**

In order to understand the diagnostic potential of the Atacama Large Millimeter/submillimeter Array (ALMA) for observations of the Sun, realistic models and corresponding artificial observations are needed. In summary, we would use a 3D model of the solar atmosphere, use it as input for a radiative transfer (RT) code, and then apply ALMA's instrumental properties to the RT results. Changing the parameters of the last step (e.g., the temporal resolution of the observation) then allows for investigating the impact of the change on the artificial observations and the extent to which initially contained information can be recovered from the artificial observations. The fundamental strength of our approach is that we can compare the processed artificial observations with the original solar input model. Such comparisons directly reveal, which aspects of the input model can be recovered or not, to what extent, and what is needed to get the most out of the data. The resulting strategies can be applied directly to real ALMA observations of the Sun and yield the potential of greatly boosting ALMA's diagnostic potential, promising high-impact results. In particular, we plan to look for observational signatures of small-scale, intermittent events in the solar atmosphere such as, e.g., various oscillations and wave modes and nano-flares, which could contribute to the still not fully explained heating of the outer layers of our Sun. Finding such small-scale events demands maximizing the capabilities of ALMA's solar observing modes. This task can be substantially supported by the numerical simulations described here but requires detailed RT calculations, which currently do not exist.  

Suitable 3D models of the solar atmosphere have been produced with the Bifrost code in Oslo and will serve as input for a new radiative transfer code that calculates how the model looks like at the wavelengths observed with ALMA. As a first step, we plan to simulate an ALMA observation run in all detail. That includes observations in two ALMA receiver bands (i.e., frequency/wavelength ranges), which each can have around 4000 spectral channels. ALMA observes at very high cadence, which was 2s in the current observing cycle but will be 1s in the cycle starting in October 2017. We plan therefore to process a time series from the 3D Bifrost model with 1s cadence over a period of 30min. In total, the radiative transfer code has to compute the emergent intensity at 2x4000 wavelengths for 1800 time steps. This task is computationally expensive and practically impossible with previously available radiative transfer codes that required substantial human interaction. The new code developed for this task has the potential to handle the challenging computations, given that it can be sped up by a factor of ideally 4, which we consider perfectly possible.  

Among other things, we plan to address the following points with the anticipated data set: 1) Tests and improvements of ALMA data calibration procedures. 2) Multi-frequency synthesis (boosting spatial resolution by sacrificing spectral resolution). 3) Novel "burst mode" (combining consecutive time steps for higher image fidelity). 

Finally, it should be noted that the exact same tools can be applied to models of other stars, too. An additional simple integration over an artificially observed stellar disk could then be used to optimize stellar observations with ALMA. Studying the solar-stellar connection with ALMA yields fundamental insights in the working of our Sun and stars in general. 



## 3. Computer resources requested  (for Types A, B, C), or expected for Type D

	Total storage required (Gbyte) 
	(only available during the duration of the preparatory access project)	
	Maximum amount of memory per core (Mbyte)	**~1GB**

## 4. Please provide the details listed below for the main simulation application  (for Types A, B, C, D)

	Name and version:  **RTam** 
	         Webpage:  **none**

	License
	If the code is open
	source please, fill out
	open source for this
	query.	

## 5. Describe the main algorithms and how they have been implemented and parallelized  
	(Maximum 250 words)

## 6. Current and target performance  (for Types A, B, C, D ; including the points below. Maximum 250 words)
	Describe the scalability of the application and performance of the application
	What is the target for scalability and performance? (i.e. what performance is needed to reach the envisaged scientific goals)

## 7. Confidentiality  (for Types A, B, C, D)
	Is any part of the project covered by confidentiality?   Yes    No
	If YES, please specify which aspect is confidential and justify:

## 8. Describe the I/O strategy regarding the parameters indicated below  (for Types A, B, C, D)
	a) Is I/O expected to be a bottleneck? (Maximum 50 words) 
	b) Implementation: I/O libraries, MPI I/O, netcdf, HDF5 or other approaches. (Maximum 50 words) Parallel HDF5
	c) Frequency and size of data output and input (Maximum 50 words) 
	d) Number of files and size of each file in a typical production run (Maximum 50 words)

## 9. Main performance bottlenecks  (for Types B, C, D. Maximum 250 words)
	

## 10. Describe possible solutions you have considered to improve the performance of the project   
	(for Types B, C, D. Maximum 250 words)
	I would go here for hybrid code (MPI + OpenMP), AVX512 and GPGPU


## 11. Describe the application enabling/optimization work that needs to be performed to achieve the target performance  (for Types B, C, D. Maximum 250 words)
	Perormance analysis + code adaptation to GPGPU/KNL 

## 12. Which computational performance limitations do you wish to solve with this project?   
	(for Types B, C, D. Maximum 250 words)
	Reduce memorry bound computations

## 13. Describe the impact of the optimization work proposed?  (for Types B, C, D. Maximum 250 words)
	This should significantly improve performance and allow us to run medium size problem on workstation size machines. 

	* Is the code widely used?
	* Would the code be used only within this original research project?
	* Would the code be used for other similar research projects with minor modifications?
	* Would the code be used in many research projects of the research field indicated in the proposal?
	* Would the modification be easy to add to the main release of the software?

## 14. Describe the request plans for work with support from PRACE experts  (for Type C, D)
	This we going to ask Marcin 

### a. Describe the level of collaboration with PRACE experts you have planned for and how much effort (person months) have you reserved for this?
	I think I can dedicate to that 4PM but falf-time and we will ask Marcin what we should put here

	Specify a rough estimate for the amount of person months this work entails :  4PM

### b. Describe the optimization work you expect to be done with the support of a PRACE expert for your project
	I think he can do OpenMP/MPI and definetly AVX512 ..I can do GPGPU prototyping

### c. Please specify the amount of PRACE experts
	person months required to support your project (1-6 PMs):  4MPs

## 15. Please describe your hardware architecture currently in use to allow an optimal selection of a feasible Tier-1 system  (for Type D only, Maximum 250 words)
	We do have Haswell but we are waiting for 2x p100 

 
