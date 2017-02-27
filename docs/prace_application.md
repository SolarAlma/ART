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

The Atacama Large Millimeter/submillimeter Array (ALMA), which is currently the worlds' largest astronomical observatory, opened up a new window on the universe. Already within the first few years of operation ALMA led to many scientific discoveries. Since December 2016, ALMA is also used for observations of our Sun. It observes the Sun at a spatial resolution, which is unprecedented in this wavelength range, and offers novel means of determining the properties of the plasma in the Sun's outer atmospheric layers. These new diagnostic capabilities, which will ultimately progress our understanding of the Sun, still need to be developed and understood further in order to fully exploit their potential. Detailed numerical simulations of the solar atmosphere and artificial observations of it play a key role in this respect. 

## 2. Scientific case of the project
Explain the scientific case for which you intend to use the code(s). Maximum 500 words.

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

 
