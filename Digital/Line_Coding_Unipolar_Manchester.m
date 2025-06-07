clear;clk;close all;
%% Initialization and Time Vector Generation
stream=randn(1,64);
bit_rate=1e4;
tb=1/bit_rate;
samples_per_rate=1000;  

num_bit=length(stream);
total_samples=num_bit * samples_per_rate;
T = num_bit / bit_rate;

Fs = samples_per_rate * bit_rate;
ts = 1 / Fs;

df=Fs/total_samples;

N = ceil(T / ts);

t = linspace(0, num_bit / bit_rate, total_samples);
if mod(N,2)==0
    f=-(Fs/2):df:(Fs/2)-df;
else
    f=-(0.5*Fs -0.5*df):df:(Fs/2)-0.5*df;
end


%% unipolar_nrz line coding  
unipolar_nrz=zeros(1,total_samples);

for i =1:num_bit
    start_index=(i-1)*samples_per_rate +1;
    end_index=i*samples_per_rate;
    unipolar_nrz(start_index:end_index) = (stream(i)>0);
end

figure;
plot(t,unipolar_nrz, 'LineWidth', 1.5);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Unipolar NRZ Signal');
xlim([0 T+5*ts]);
ylim([-0.2 1.2]);             
grid on;    


%% Manchester line coding 
Manchester=zeros(1,total_samples);
for i =1:num_bit
    start_index=(i-1)*samples_per_rate +1;
    end_index=i*samples_per_rate;
    mid_index = start_index + floor(samples_per_rate/2) - 1;
    if(stream(i)>0)
        Manchester(start_index:mid_index)=1;
        Manchester(mid_index+1:end_index)=-1;
    else
        Manchester(start_index:mid_index)=-1;
        Manchester(mid_index+1:end_index)=1;
    end
end    
figure;
plot(t,Manchester, 'LineWidth', 1.5);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Manchester Signal');
xlim([0 T+5*ts]);
ylim([-1.2 1.2]);             
grid on;    



%% FFT and PSD for Unipolar NRZ
M_unipolar = fftshift(fft(unipolar_nrz))/N;

unipolar_psd = abs(M_unipolar).^2;  


figure;
plot(f*tb, abs(M_unipolar), 'LineWidth', 1.5);
xlabel('Nomarlized Frequency (Hz)');
ylabel('Magnitude');
title('FFT Magnitude of Unipolar NRZ Signal');
xlim([-3 3]);
grid on;


figure;
plot(f*tb, unipolar_psd, 'LineWidth', 1.5);
xlabel('Nomarlized Frequency (Hz)');
ylabel('PSD');
title('PSD of Unipolar NRZ Signal');
xlim([-1 1]);
grid on;

%% FFT and PSD for Manchester

M_Manchester = fftshift(fft(Manchester))/N;

Manchester_psd = abs(M_Manchester).^2;  


figure;
plot(f*tb, abs(M_Manchester), 'LineWidth', 1.5);
xlabel('Nomarlized Frequency (Hz)');
ylabel('Magnitude');
title('FFT Magnitude of Manchester Signal');
xlim([-5 5])
grid on;


figure;
plot(f*tb, Manchester_psd, 'LineWidth', 1.5);
xlabel('Nomarlized Frequency (Hz)');
ylabel('PSD');
title('PSD of Manchester Signal');
xlim([-5 5])
grid on;




