## Test environments
* local OS X install, R 3.6.0
* Ubuntu 16.04.6 LTS (on travis-ci), R 4.0.2 
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE for win-builder:

* checking CRAN incoming feasibility ... NOTE
Possibly mis-spelled words in DESCRIPTION: 
  Deconvolution (3:19)
  Kang (12:152)
  al (12:160)
  et (12:157)
  
Deconvolution is NOT a mis-spelled word. Kang et al. is the paper that describes the method.

## Downstream dependencies
There are currently no downstream dependencies for this package

## Resubmission
This is a resubmission. In this version I have:

* Added citation information for the package.

* Added a few more examples in the Rd-files and enabled automatic tesing using usethis::use_testthat().

* Added option to suppress printing message to console, i.e. setting verbose = FALSE.

* Added Authors@R field.

* Reset par() option to default after plotting.

* Ensured examples do not use more than 2 cores.

## Resubmission
* Added all authors, contributors and copyright holders in the Authors@R field with the appropriate roles.

## Resubmission 2020-02-05
* Added verbose parameter to prevent from writing to console.

* Makoto Matsumoto and Takuji Nishimura are NOT part of the R package or method development. I used the Mersenne Twister random number generator developed by Makoto Matsumoto and Takuji Nishimura. The Mersenne Twister random number generator code is in cokus.cpp and I quote the writer of the cokus.cpp file: "This version is a recode by Shawn Cokus (Cokus@math.washington.edu) on March 8, 1998 of a version by Takuji Nishimura (who had suggestions from Topher Cooper and Marc Rieffel in July-August 1997)."





