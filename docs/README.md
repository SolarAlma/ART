prace_application.md


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




