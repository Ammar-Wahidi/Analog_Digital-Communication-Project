# Analog and Digital Communication Project

## Overview
This repository contains the code and documentation for the **Analog & Digital Communication Project** submitted for the course **Fundamental of Communication Systems (ECE 252)**. The project is divided into two main parts: **Analog Communication** and **Digital Communication**, focusing on signal processing, modulation, and demodulation techniques using Matlab.

The project was completed by:
- Amr Ahmed Wahidi 
- Ammar Ahmed Wahidi 
- Abdelrahman Essam Fahmy 
- Mohamed Yehia Saeed 
- Sobhey Mohamed Osman 

**Submitted to**: Dr. Mazen Erfan, Dr. Alaa Fathy, Eng. Ahmed Al-Sayed, Eng. Ahmed Khaled  
**Submission Date**: June 5th, 2025

## Project Structure
The project is organized into two primary sections, each addressing different aspects of communication systems:

### 1. Analog Communication
This section focuses on the analysis, modulation, and transmission of analog signals using techniques such as:
- Signal generation and time-domain plotting.
- Fourier transform derivation and FFT computation.
- Bandwidth estimation based on power spectral density.
- Low-pass filtering (LPF) with bandwidths of 1 Hz and 0.3 Hz.
- Modulation schemes: DSB-SC (Double-Sideband Suppressed Carrier), SSB (Single Sideband, USB), and FDM (Frequency Division Multiplexing).
- Coherent detection for signal recovery.

Key tasks include:
- Plotting the function \( x(t) \) and its Fourier transform.
- Filtering and reconstructing signals.
- Implementing DSB-SC and SSB modulation with a 2 Hz guard band.
- Visualizing FDM signals and performing coherent demodulation.

### 2. Digital Communication
This section explores digital signal processing and modulation techniques, including:
- **Line Coding**: Implementation and comparison of Unipolar NRZ and Manchester coding for a 64-bit random stream.
- **Binary Phase-Shift Keying (BPSK)**: Modulation and coherent demodulation of a random bit stream, including analysis of phase offsets (30°, 60°, 90°) and bit error rate (BER).

Key tasks include:
- Generating and visualizing coded signals in the time domain.
- Computing and plotting FFT and Power Spectral Density (PSD) for line codes.
- Implementing a BPSK transmitter and receiver, analyzing the spectrum, and evaluating performance under different phase conditions.

## Repository Contents
- **/analog**: Contains Octave scripts for analog communication tasks, including signal plotting, Fourier transforms, filtering, and modulation (DSB-SC, SSB, FDM).
- **/digital**: Contains Octave scripts for digital communication tasks, including line coding (Unipolar NRZ, Manchester) and BPSK modulation/demodulation.
- **/report**: The project report (`Communication Project Report.pdf`) detailing the methodology, code, results, and analysis.
- **README.md**: This file, providing an overview of the project and instructions.

## Prerequisites
To run the scripts in this repository, you need:
- **MATLAB** (or Octave) installed on your system.
- Basic understanding of signal processing and communication systems.

## How to Run
1. Clone the repository:
   ```bash
   git clone https://github.com/Ammar-Wahidi/Analog_Digital-Communication-Project.git
   ```
2. Navigate to the `analog` or `digital` directory to access the respective scripts.
3. Open Octave and run the scripts (e.g., `script_name.m`) to execute the tasks and generate plots.
4. Refer to the project report (`/report/Communication Project Report.pdf`) for detailed explanations of each task and results.

## Key Features
- **Analog Communication**:
  - Time-domain and frequency-domain analysis of signals.
  - Implementation of DSB-SC, SSB (USB), and FDM modulation schemes.
  - Bandwidth estimation and low-pass filtering.
  - Coherent detection for signal recovery.
- **Digital Communication**:
  - Comparison of Unipolar NRZ and Manchester line coding in time and frequency domains.
  - BPSK modulation and demodulation with phase offset analysis.
  - Visualization of signals, FFT, PSD, and bit error rates.

## Results
- **Analog Communication**: Successfully plotted and analyzed signals, derived Fourier transforms, estimated bandwidth, and implemented modulation/demodulation schemes. The FDM scheme effectively combined DSB-SC and SSB signals with a 2 Hz guard band.
- **Digital Communication**: Demonstrated the differences between Unipolar NRZ and Manchester coding, highlighting Manchester’s superior clock recovery and DC balance. BPSK modulation showed robust performance, with zero bit errors at 30° and 60° phase offsets, but complete mismatch at 90° due to orthogonality.

## References
- *Introduction to Analog and Digital Communications* by Simon Haykin and Michael Moher.

