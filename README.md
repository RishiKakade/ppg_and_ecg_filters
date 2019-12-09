# ppg_and_ecg_filters
A collection of subroutines and functions for removing noise and processing ECG and PPG signals I found useful.

Some background, my ECG sensor had a sampling frequency of 250hz and my [PPG](https://www.nonin.com/products/oem3/) had one of 75hz 
### Transfer Functions

##### Chebyshev Lowpass
  - Order 2

  
##### Butterworth Lowpass
  - Order 5
  - Baseline removal (sub 0.5-0.67Hz)

##### Butterworth Highpass
  - Order 2
  - HF removal (150Hz)

##### Butterworth Bandpass
  - Order 2
  - Reasonable passband is 0.5-150Hz. [Here's why](https://www.researchgate.net/publication/303155606_A_survey_of_noise_removal_techniques_for_ecg_signals)

##### Butterworth Bandstop
  - Order 2
  - In North America, the electrical grid operates at 60Hz

### Other Processing Tools
##### QRS Peak Detection
  - My method of finding R peaks is pretty much returns indicies with ECG measurments of >750uV. Theres not really any sense in applying complicated peak detection algs if its a clean signal

##### Saturation Correction
  - My ECG had an upper limit of 50k uV

##### Moving Median Filter
  - This one is more useful for PPG signals. You could use a bandpass+notch but its not realy necessary if you're just trying to compute heart rate or something. The R peak detection algorithm works here too, just with a different bound 

##### BPM and HR Variance
  - BPM is just the number of R peaks per minute. You could use Q or S peaks too I guess, but R peaks are more prominent thus likely give a more accurate heart rate
  - Variance just uses a built in numpy function. This is more useful on time frames where you can approximate as roughly linear. Variance doesn't tell you anything if a person jumps out of bed and goes for crossfit!
