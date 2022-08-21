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
1. **Main file** (SSS_DeepLearning.m): This is the main file to execute.
2. **Inputs:**
  **Sound Data**
    Two .WAV files (equipment sound collected on the job site using a microphone)
      _ori: means original sound
      _den: means denoised sound

## Methodology
![image](https://user-images.githubusercontent.com/73087167/185807793-cc696857-14ce-40c0-8682-3ca6ca62cc0f.png)


## Results (Hard TFM)
![image](https://user-images.githubusercontent.com/73087167/185807809-dd383387-a636-4b70-8f94-c43b92f91833.png)

## Results (Soft TFM)
![image](https://user-images.githubusercontent.com/73087167/185807818-a4dbf015-ef4f-48c3-bbcc-eb753f3cb6f2.png)

## Source Separation Evaluation
![image](https://user-images.githubusercontent.com/73087167/185807835-7300e7c0-87ae-4d67-a15b-0183756c45e9.png)
