This folder contains frequencies used for Anavar inference of DFEs. Because Anavar cannot use a truncated SFS, we have to extend the frequency bins from 1-98%, into 50 bins.
We simulated data to correct the SFS due to different levels of sequencing errors.

* anavar/ run without error calibration
* anavar_calibrated_err0.01/ calibrated based on simulations with sequencing error of 0.01
* anavar_calibrated_err0.0033/ calibrated based on simulations with sequencing error of 0.0033
* anavar_calibrated_err0.001/ calibrated based on simulations with sequencing error of 0.001
* Check the README.txt files in anavar/ and anavar_calibrated_err0.001/. The other two scenarios are similar.
