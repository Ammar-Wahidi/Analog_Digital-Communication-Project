clc;clear; close all ;
%% Degien Terms
df= 0.01;
fs = 100;
ts = 1/fs;
T = 1/df;
N = round(T*fs);
t = -T/2: ts: T/2-ts;
f = (-fs/2):df:(fs/2-df);

%% Plot the function x(t) shown on Octave 
 x = zeros(size(t));
 x(t>=-4 & t<0) = t(t>=-4 & t<0) + 5;
 x(t>=0 & t<=4) = -t(t>=0 & t<=4) + 5;
 figure(1);
 plot (t,x);
 xlabel('t');
 ylabel('x(t)');
 title('message x(t) in time domain');
 xlim([-5 5]);
 grid on;                %question 1
%% Derive an analytical expression for its Fourier transform
G = 8*sinc(8*f)+16*(sinc(4*f).^2); %analytical
%% FFT of the signal 
XF = fftshift(fft(x))*ts;
XF_mag = abs(XF);          % function is non-periodic so multiply with ts
figure;
plot(f,XF_mag, 'b');        
hold on;
plot (f,abs(G), 'r');
xlabel('frequency');
ylabel('magnitude')
title('fourier transform comparison');
legend('FFT result','Analytical');
xlim([-3 3]);
grid on;
%% Estimate the BW defined as the frequency band after which the power spectrum of the signal drops to 5% of its maximum value.
power_m = abs(XF.^2); 
P_mag_max = max(power_m); 
threshold = 0.05 * P_mag_max; 
index1 = find(power_m >= threshold);
BW1_estimate = max(abs(f(index1)))
%%  LPF Filtering with BW=1Hz
BW = 1;
H_lpf1 = abs(f) < BW; % Renamed H to avoid conflict
x_rec = real(ifft(ifftshift(H_lpf1.*XF)/ts));
figure;
plot(t,x_rec);
hold on;
plot(t,x);
legend('after LPF','input signal');
title('x(t) reconstruction after LPF (BW=1Hz)');
xlabel('time');
ylabel('Amplitude');
xlim([-15 15]);
ylim([-1 6]);
grid on;
%% LPF Filtering with BW=0.3Hz
BW_03 = 0.3; % Renamed BW to avoid conflict
H_lpf2 = abs(f) < BW_03; % Renamed H to avoid conflict
x_rec2 = real(ifft(ifftshift(H_lpf2.*XF)/ts));
figure;
plot(t,x_rec2,'b');
hold on;
plot(t,x,'r');
legend('after LPF=0.3Hz', 'input signal');
title('x(t) reconstruction after LPF (BW=0.3Hz)');
xlabel('time');
ylabel('X(t)');
xlim([-15 15]);
ylim([-1 6]);
grid on;
%% Analysis of m(t) = cos(2π·0.5·t) for 0 < t < 4
% Step 1
m = zeros(size(t));
m(t>0 & t<4) = cos(2*pi*0.5*t(t>0 & t<4));
figure;
plot(t,m);
title('Message m(t) in time domain.');
xlabel('time');
ylabel('m(t)');
xlim([-10 10]);
ylim([-1.2 1.2]);
grid on;
% Step 2
GF= 2*sinc(4*(f-0.5))+ 2*sinc(4*(f + 0.5));
% Step 3
M = fftshift(fft(m))*ts;
figure;
plot(f,abs(M),'b');
hold on;
plot(f,abs(GF),'r');
title('Fourier transform comparison for m(t)');
legend('FFT result','Analytical');
xlabel('frequency');
ylabel('Magnitude');
xlim([-12 12]);
grid on;
% Step 4
power_m = abs(M.^2); 
P_mag_max = max(power_m); 
threshold = 0.05 * P_mag_max; 
index2 = find(power_m >= threshold);
BW2_estimate = max(abs(f(index2)))

%% FDM Modulation Scheme
% Modulation of message x(t) using DSB-SC modulation
fc1 = 20;
carrier1 = cos(2*pi*fc1*t);
s1= x_rec.*carrier1; % Using x_rec as presumably the filtered version is intended for modulation
S1 = fftshift(fft(s1))*ts;
figure;
plot(f,abs(S1));
xlabel('f (Hz)');
ylabel('|S1(f)|');
title('DSB-SC modulation of x_{rec}(t)');
xlim([-30 30]);
grid on;

% Modulation of message m(t) using SSB modulation
fc2 = 23;               %USB is used
carrier2 = cos(2*pi*fc2*t);
s2 = m.*carrier2;
S2 = fftshift(fft(s2))*ts;
H_ssb = zeros(size(f)); % Renamed H to avoid conflict
% For USB, passband is [fc2, fc2 + BW_of_m] and [-fc2 - BW_of_m, -fc2]
% Assuming BW_of_m is around 1 Hz for m(t) = cos(pi*t) for 0<t<4
BW_m_approx = 1; % Approximate bandwidth of m(t)
H_ssb(f>=fc2 & f<=(fc2+BW_m_approx))=1; % Positive frequencies for USB
H_ssb(f<=-fc2 & f>=(-fc2-BW_m_approx))=1; % Negative frequencies for USB

S2SSB = S2 .*H_ssb;
s2ssb = real(ifft(ifftshift(S2SSB)/ts)); % Use real for ifft of SSB signal
figure;
plot(f,abs(S2SSB));
title('SSB-USB modulation of m(t)');
xlabel('f (Hz)');
ylabel('|S2_{SSB}(f)|');
xlim([-30 30]);

grid on;

stotal = s1 + s2ssb; % FDM: s1 (DSB-SC) + s2ssb (SSB)
figure;
plot(t,stotal);
xlabel('t (seconds)');
ylabel('s_{total}(t)');
title('FDM signal s_{total}(t) in time domain');
xlim([-25 25]);
grid on;

%% FDM Signal Visualization
Stotal = fftshift(fft(stotal))*ts;
figure;
plot(f,abs(Stotal));    %two spectrums on the same graph
xlabel('f (Hz)');
ylabel('|S_{total}(f)|');
title('FDM signal S_{total}(f) in frequency domain');
xlim([-35 35]);
grid on;

%% Coherent Detector
% for message x(t):
x_before_LPF = s1 .*carrier1; % s1 already is x_rec . carrier1
X_before_LPF = fftshift(fft(x_before_LPF))*ts;
figure;
plot(f,abs(X_before_LPF));
title('Spectrum of x_{rec}(t) after mixing with carrier (before LPF)');
xlabel('f (Hz)');
ylabel('Magnitude');
grid on;

H_demod_x = abs(f)<BW; % Using BW=1 for LPF for x(t) (consistent with x_rec)       
x_received = 2*real(ifft(ifftshift(H_demod_x.*X_before_LPF)/ts)); % Multiply by 2 for DSB-SC demod
figure;
plot(t,x,'r'); % Plot original x for comparison
hold on;
plot(t,x_received,'b');
xlabel('time (seconds)');
ylabel('signal x(t)');
title('x(t): original & received (after DSB-SC demod)');
legend('original signal','received signal');
xlim([-6 6]);
grid on;

% for message m(t)
% For SSB, coherent detection involves multiplying by carrier and LPF
m_before_LPF = s2ssb .*carrier2; 
M_before_LPF = fftshift(fft(m_before_LPF))*ts;
figure;
plot(f,abs(M_before_LPF));
title('Spectrum of m(t) after mixing with carrier (SSB demod, before LPF)');
xlabel('f (Hz)');
ylabel('Magnitude');
grid on;

% Bandwidth for LPF for m(t) should be related to BW_m_approx
H_demod_m = abs(f)<BW_m_approx; 
m_received = 2*real(ifft(ifftshift(H_demod_m.*M_before_LPF)/ts)); % Multiply by 2 for SSB demod with cos
figure;
plot(t,m,'r');
hold on;
plot(t,m_received,'b');
xlabel('time (seconds)');
ylabel('signal m(t)');
title('m(t): original & received (after SSB demod)');
legend('original signal','received signal');
xlim([-15 15]);
ylim([-1.1 1.1]);
grid on;