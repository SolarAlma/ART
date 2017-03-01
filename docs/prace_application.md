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

 
===
===
# SUBMITTED VERSION (to be formatted) 
===

Export to :
Excel   PDF  
PRACE banner : Partnership for Advanced Computing in Europe Front page  > Proposal 2010PA3776 edition Sign out
Proposal n° 2010PA3776 – Application form

Your proposal has successfully been created.

If you make modifications on this page and do not press the save button at the bottom of the page then your changes will not be saved. Pressing save does not submit your proposal. To submit your proposal, you must follow the link on the proposals summary page. After submission, you may make changes to your proposal. You will need to un-submit your proposal, make the necessary changes, save those changes and submit again your proposal. Please note that proposals not submitted will not be assessed. All information must be given in English.

Please note that mandatory fields are indicated by a red square (required).
General information
Type of proposal:	required	
	A – Code scalability testing.
Scalability testing to obtain scalability plots which can be used as supporting information
when applying to future PRACE project.
Please fill out parts 1-8 of the application form
	B – Code development and optimization by applicant (without PRACE support).
Please fill out parts 1-13 of the application form.
	C – Code development with support from PRACE experts.
Please fill out parts 1-14 of the application form
	D – Code development, starting on Tier-1 system with support from PRACE experts.
Please fill out parts 1-15 of the application form
 
Start date:		
01.05.2017
 
The start date should be either one and half or two months
after the announced date for evaluation cut-off.
 
Is this proposal a :	 	
 Continuation of Type D project (only for Type A)
ID of the previous proposal :  
2010PAxxxx
Project name:	required	
Radiative Transfer Forward Modelling of Solar Observations with ALMA
Research field:	required	
Project leader (personal data and contact)
Gender:	required	 Female  Male
Title:	required	
Dr.
 e.g. Mr., Mrs, Miss, Dr, Prof.
First name:	required	
Sven
Last name:	required	
Wedemeyer
Initials:		
SW
Date of birth:	required	
20061975
 ddmmyyyy, e.g. 01021970
Nationality:	required	
e-mail :	required	svenwe@astro.uio.no
Phone number (direct):	required	
+4722856520
 Prefix '+' followed by country code, e.g. +33 for France.
Project leader (organisation and job title)
Job title:	required	
Researcher/Project leader
Website:		
https://folk.uio.no/svenwe/
Organisation name:	required	
University of Oslo
Department:	required	
Institute of Theoretical Astrophysics
Group:		
Solar Physics
Address:	required	
Postboks 1029 Blindern

Postal code:	required	
0315
City:	required	
Oslo
Country:	required	
Organisation with a research activity: required   Yes     No
Employment contract of the project leader is valid at least 3 months after the end of the allocation period: required   Yes     No
For commercial companies,
Is the head office of the organisation in Europe?    Yes     No
% of R&D activity in Europe as compared to total R&D activity :  
 Contact person for all correspondence
Please give your professional e-mail address.
E-mail addresses such as Gmail and Hotmail are not accepted.

Name:	required	
Mikolaj Szydlarski
E-mail:	required	
mikolaj.szydlarski@astro.uio.no
 Requested computer systems required 
You may choose more than one computer system.

Computing center	Machine
BSC
Barcelona Supercomputing Center 	 MareNostrum	
BSC
Barcelona Supercomputing Center 	 MareNostrum Hybrid Nodes	
CINECA
 	 Marconi – Broadwell	
CINECA
 	 Marconi – KNL	
CSCS
Swiss National Supercomputing Centre 	 Piz Daint	
Gauss/HLRS
High Performance Computing Center Stuttgart (HLRS) 	 Hazel Hen	
Gauss/JSC
Juelich 	 Juqueen	
Gauss/LRZ
Leibniz-Rechenzentrum 	 SuperMUC	
GENCI/CEA
Commissariat à l'Énergie Atomique 	 Curie Thin Nodes (TN)	
   - or -
Tier-1 system (for PA Type D)
A Tier-1 system will be selected based
on your given information		
Specific Tier-1 system (for PA Type D)
Name of the preferred Tier-1 system	
Abel cluster at UiO

PLEASE NOTE that there is NO prototype system available for access at this cut-off date.

1. Summary of the project  (for Types A, B, C, D) 
To be published in the PRACE website. Maximum 500 words.


The Atacama Large Millimeter/submillimeter Array (ALMA), which is currently the world's largest astronomical observatory, opened up a new window on the universe. The interferometric array is located on the Chajnantor plateau in the Chilean Andes at an altitude of 5000 m and consists of 66 antennas, most of them with a diameter of 12 m. By combining the antennae, they act like a giant telescope with baselines of up to 16 km. Already within the first few years of operation, ALMA led to many scientific discoveries. Since December 2016, ALMA is also used for observations of our Sun. It observes the Sun at a spatial resolution, which is unprecedented in this wavelength range, and offers novel means of determining the properties of the plasma in the Sun's outer atmospheric layers. Due to the properties of the solar radiation at millimeter wavelengths, ALMA serves as a linear thermometer, mapping narrow layers at different heights. It can measure the thermal structure and dynamics of the solar atmosphere and thus sources and sinks of atmospheric heating.

Among other expected scientific results, ALMA promises significant steps towards understanding the intricate dynamics and physical processes that, in combination, might yield the solution of the coronal heating problem - a long standing fundamental question in modern astrophysics. However, ALMA's novel diagnostic capabilities, which will ultimately progress our understanding of the Sun, still need to be developed and understood further in order to fully exploit the instrument’s potential. Detailed numerical simulations of the solar atmosphere and artificial observations of the Sun play a key role in this respect.

Such artificial observations of the Sun will be produced as part of the SolarALMA project at the University of Oslo, which is funded with a Consolidator Grant by the European Research Council (ERC), in collaboration with Dr. de la Cruz Rodriguez from the University of Stockholm. The overall aim of the SolarALMA project is to utilize the first observations of the Sun with ALMA and to develop the required observing and analysis techniques. An important step in this endeavour is the development of realistic numerical models of the solar atmosphere, which can be used to test how to optimally set up solar observations, and how to analyse and interpret them. While 3D numerical models are routinely produced on high-performance computers, the codes available for producing the corresponding (artificial) ALMA observations of such models did not perform well enough so far. We have developed a new code that solves the radiative transfer equation for a 3D numerical model and thus reveals how the modeled part of the Sun would look like through the eyes of ALMA at (sub)millimeter wavelengths. The new code is in an advanced stage but, still needs to be optimized in order to provide the basis for essential studies, which will result in optimizing ALMA's scientific impact for observing and understanding the Sun.
2. Scientific case of the project  (for Types A, B, C, D) 
Explain the scientific case for which you intend to use the code(s). Maximum 500 words.


In order to understand the diagnostic potential of the Atacama Large Millimeter/submillimeter Array (ALMA) for observations of the Sun, realistic models and corresponding artificial observations are needed. In summary, we would use a 3D numerical model of the solar atmosphere as input for a radiative transfer (RT) analysis, which outputs corresponding images at mm wavelengths, and then apply ALMA's instrumental effects (e.g., limited resolution, noise) to these images. Changing the instrumental parameters in the last step  (e.g., the temporal resolution) then allows for investigating the impact of the change on the processed images and the extent to which initially contained information (e.g., temperature variations on small spatial scales) can be recovered from these "simulated observations". The fundamental strength of our approach is that we can compare the simulated observations with the original solar atmosphere  model. Such comparisons directly reveal, which aspects of the input model can be recovered and what is needed to get the most out of the data. The resulting strategies can be applied directly to real ALMA observations of the Sun, potentially greatly boosting ALMA’s diagnostic capabilities. In particular, we plan to look for observational signatures of small-scale, intermittent events in the solar atmosphere such as, e.g., oscillations, waves, and nano-flares, which could contribute to the still not fully explained heating of the outer solar layers. Finding such small-scale events demands maximizing the capabilities of ALMA's solar observing modes. This task can be substantially supported by the numerical simulations described here, but requires detailed RT calculations, which currently do not exist. 

Suitable 3D models of the solar atmosphere have been produced with the Bifrost code in Oslo and will serve as input for a new radiative transfer code that calculates how the model looks like at the wavelengths observed with ALMA. As a first step, we plan to simulate an ALMA observation run in detail. That includes observations in two ALMA receiver bands (i.e., wavelength ranges) with 4000 spectral channels each. ALMA observes at very high cadence, which will be 1s in the cycle starting in October 2017. We plan therefore to process a time series from the 3D Bifrost model with 1s cadence over a period of 30min. In total, the radiative transfer code has to compute the emergent intensity at 2x4000 wavelengths for 1800 time steps. This task is computationally expensive and practically impossible with previously available radiative transfer codes that required substantial human interaction. The new code developed for this task has the potential to handle the challenging computations, given that it can be sped up by a factor of ideally 4, which we consider perfectly possible. 

Among other things, we plan to address the following points with the anticipated data set: 1) Tests and improvements of ALMA data calibration procedures. 2) Multi-frequency synthesis (boosting spatial resolution by sacrificing spectral resolution). 3) Novel "burst mode" (combining consecutive time steps for higher image fidelity).

Finally, it should be noted that the same tools can be applied to models of other stars, too.
3. Computer resources requested  (for Types A, B, C), or expected for Type D

Total storage required (Gbyte)
(only available during the duration of
the preparatory access project)	
1000
Maximum amount of memory per core (Mbyte)	
1000
4. Please provide the details listed below for the main simulation application  (for Types A, B, C, D)

Name and version	
RTma/Nyx 
(version: alpha) 
Webpage
or other reference	
none yet 
License
If the code is open
source please, fill out
open source for this
query.	
It is currently an alpha version without license assigned but we plan to release it under GPL or MIT license.
5. Describe the main algorithms and how they have been implemented and parallelized  
(for Types A, B, C, D ; Maximum 250 words)


The code is written in C++ with some features from C++11 standard with one additional module written in Fortran77.

At the core of computation is a nonlinear solver for radiative transfer (RT), which requires the solution of the (polarized) radiative transfer problem. This problem involves the solution of a set of coupled first-order linear differential equations, one for each Stokes parameter.

Given that each column in the atmosphere model is a completely independent calculation, the problem is very suitable for parallel computing. A very limited amount of communication is necessary between different processes. The code has been parallelized using MPI. The total number of columns to be calculated (tasks) is divided by the number of processes, so that the amount of tasks each process has is about the same and is known in the beginning of the execution. Each process starts working through its task list until it is finished, and then writes the output to the disk and waits for the other processes to complete.

6. Current and target performance  (for Types A, B, C, D ; including the points below. Maximum 250 words)

Describe the scalability of the application and performance of the application
What is the target for scalability and performance? (i.e. what performance is needed to reach the envisaged scientific goals)

The code is naively parallel but it is very computationally demanding. It scales with ease up to 512 cores but it has not been tested yet  with a higher number of CPUs. We believe that better performance can be achieved with a hybrid approach (OpenMP/MPI) or by using accelerators like Intel Xeon Phi or GPGPU.

We would need a speed up by a factor of ideally 4 to be able generate the artificial observations for ALMA as described above (item 2).
7. Confidentiality  (for Types A, B, C, D)

Is any part of the project covered by confidentiality?   Yes    No

If YES, please specify which aspect is confidential and justify:


8. Describe the I/O strategy regarding the parameters indicated below  (for Types A, B, C, D)

a) Is I/O expected to be a bottleneck? (Maximum 50 words)


Although the input and output is very volumetric, the time spend on I/O operations is negligible when compared to the computation time.
b) Implementation: I/O libraries, MPI I/O, netcdf, HDF5 or other approaches. (Maximum 50 words)


We use the parallel HDF5 library for I/O.
c) Frequency and size of data output and input (Maximum 50 words)


Input read only at the beginning (~25GB for snapshot),  output saved only at the end of computation (~20GB for each snapshot)
d) Number of files and size of each file in a typical production run (Maximum 50 words)


1 hdf5 file of 1 model snapshot (~25GB) for input and 1 hdf5 for output (~20GB). 

The code can process time series of model snapshots sequentially. The large production run mentioned under item 2 would thus result in 1800 I/O files (processed one at a time).  
9. Main performance bottlenecks  (for Types B, C, D. Maximum 250 words)


Currently the slowest part is the EOS (Equation of State) module which strongly depends on expensive calculation of integrals of exponential functions.
10. Describe possible solutions you have considered to improve the performance of the project   
(for Types B, C, D. Maximum 250 words)


The major improvement would require changes in the algorithm. We believe that using interpolation tables for the equation of state instead of the explicit calculations of integrals of exponential functions would speed up the code substantially keeping reasonable precision of the solution.

The code has not been optimized so far for better vectorization, thus in depth performance analysis is needed with fine tuning (loop by loop) for better SIMD usage. We also believe that our code should be easy to accelerate on Xeon Phi (KNL) or/and GPGPU.  

Adding OpenMP should also reduce memory overhead and would improve scalability.

11. Describe the application enabling/optimization work that needs to be performed to achieve the target performance  (for Types B, C, D. Maximum 250 words)


1) Vectorization guidance (e.g., via OpenMP 4.0, compiler pragmas, or intrinsics)
2) Investigation of code suitability and code adaptation to GPGPU/KNL
3) Introducing OpenMP and multi-threading to reduce overheads on many-core, shared memory architectures

 

12. Which computational performance limitations do you wish to solve with this project?   
(for Types B, C, D. Maximum 250 words)


The goal is to improve the performance of the EOS module and other computationally intense parts of the code, in particular

 - Reduce number of memory bound computations
 - Improve vectorization (increase FLOPs)
 - Validate whether it would be possible to port this code to GPGPU
 - Provide performance estimate for/on Intel Xeon Phi KNL

13. Describe the impact of the optimization work proposed?  (for Types B, C, D. Maximum 250 words)

Is the code widely used?
Would the code be used only within this original research project?
Would the code be used for other similar research projects with minor modifications?
Would the code be used in many research projects of the research field indicated in the proposal?
Would the modification be easy to add to the main release of the software?

The code is not widely used yet, but many applications are possible and foreseen, not only for modelling of the Sun but also other stars. 

The code would be an integral part of the SolarALMA project. It is explicitly mentioned as a work-package in the ERC project plan. 

We intend to make the code available for a large user base, in particular for the Solar Simulations for the Atacama Large Millimeter Observatory Network (SSALMON), which is hosted and organized by our group. It currently involves 86 scientists from around the world.  

SSALMON www page: http://www.ssalmon.uio.no/
14. Describe the request plans for work with support from PRACE experts  (for Type C, D)

14.a. Describe the level of collaboration with PRACE experts you have planned for and how much effort (person months) have you reserved for this?


Our plan is to dedicate one member of our project to work on the algorithmic part and generating input data and validating the output.  
Specify a rough estimate for the amount of person months this work entails :  
6
14.b. Describe the optimization work you expect to be done with the support of a PRACE expert for your project


While we have resources to work internally on improving the algorithm and introducing OpenMP, we would like a PRACE expert to help us with performance analysis, vectorization and optimization for our local Tier-1 system (Abel cluster at UiO) as an intermediate stage towards larger HPC systems. We also would like to work in sync on porting our code to GPGPU/KNL.


14.c. Please specify the amount of PRACE experts
person months required to support your project (1-6 PMs):   
4

15. Please describe your hardware architecture currently in use to allow an optimal selection of a feasible Tier-1 system  (for Type D only, Maximum 250 words)


 We regularly run our codes on a small resource our group has, with multi-core Intel Xeon, 64GB/compute node. Some of the nodes are equipped with Infiniband. Access to, and help with a similar, larger Tier-1 system with some modern GPUs (e.g., Nvidia p100) and Xeon KNL, where we could test the adapted/optimized the code, is required for us to achieve the scientific project goals.
 

  required I certify that I have read, understand, accept and comply with the terms and
conditions of PRACE Preparatory access – Call for proposals available at
http://prace-ri.eu/PRACE-Preparatory-Access
Those terms include the ones reproduced hereinafter for the sake of clarity:

The users commit to:

Provide to PRACE within the period established in the guide for applicants a final report, using the proper PRACE template, with the results obtained through the access to the PRACE Research Infrastructure, as well as a qualitative feedback on the use of the resources.
Acknowledge the role of the HPC Centre and PRACE in all publications which include the results above mentioned. Users shall use the following (or equivalent) wording in such acknowledgement in all such papers and other publications: 

« We acknowledge PRACE for awarding us access to resource [machine name] based in [country] at [site] » 

Where technical support has been received the following additional text should also be used: 

« The support of [name of person/people] from [organisation name], [country] to the technical work is gratefully acknowledged. » 

Allow PRACE to publish the mentioned report as of one year from the termination of the allocation period.
Commit to collaborate with PRACE, upon its request, in the preparation of dissemination material.
The applicant commits to not use the project results for military purposes.
   SAVE   
v1.6.1 © Cines.fr
Contact us here for queries regarding how to fill in the online form
PRACE website



