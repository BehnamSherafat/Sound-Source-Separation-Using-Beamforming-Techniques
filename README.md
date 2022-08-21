# Sound-Source-Separation-Using-Beamforming-Techniques
Please cite this paper if you used the codes or dataset:
Sherafat, B., Rashidi, A. and Asgari, S., 2020, December. Comparison of different beamforming-based approaches for sound source separation of multiple heavy equipment at construction job sites. In 2020 Winter Simulation Conference (WSC) (pp. 2435-2446). IEEE.
## Keywords
1. Beamforming Techniques
  Time-delay Beamformer
  Sub-band Phased Shift Beamformer
  Time-delay Linear Constraint Minimum Variance (LCMV) Beamformer
  Frost Beamformer
  Generalized Side-lobe Canceler (GSC) Beamformer
  Wide-band Minimum-Variance Distortionless-Response (MVDR) Beamformer
2. Microphone arrays
3. Direction of Arrival (DOA)
4. Off-the-shelf microphone
5. Conventional or adaptive beamformer
6. Narrowband or wideband beamformer
7. Time-domain or frequency-domain beamformer
8. Array gain

## Introduction
Six beamforming-based approaches for construction equipment sound source separation are implemented and evaluated using real construction job site data. The results show that Frost beamformer and time-delay Linear Constraint Minimum Variance (LCMV) generate outputs with array gains of more than 4.0, which are more reliable for audio-based equipment activity recognition.

## Programming Language Used
MATLAB

## Guide to folders and files
Codes are explained as follows:
1. **Main file** (EquipmentBeamformers.m): This is the main file to execute.
2. **Inputs:**
  Three types of equipment and each equipment has Two .WAV files (equipment sound collected on the job site using a microphone array)
  1. BIGD9GBulldozer (_den: means denoised sound)
  2. Excavator (_den: means denoised sound)
  3. Jackhammer (_den: means denoised sound)
      

## Delay-and-sum beamforming architecture
![image](https://user-images.githubusercontent.com/73087167/185808285-706838dd-1e34-4fb6-8f4e-d4e3aee49cf2.png)

## The filter-sum beamformer structure for Frost’s algorithm (Greensted 2010)
![image](https://user-images.githubusercontent.com/73087167/185808302-b9e70677-d568-414a-8885-3615a0c45d8d.png)

## Data Collection Setup
![image](https://user-images.githubusercontent.com/73087167/185808313-be3cf812-f015-41b7-a13d-8a0486b40e48.png)

## Comparison of different beamformers
![image](https://user-images.githubusercontent.com/73087167/185808336-84e29aec-bcf9-4850-8440-86c554ea09df.png)





