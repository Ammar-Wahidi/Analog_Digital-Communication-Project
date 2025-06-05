%% Generate a random bitstream
%stream = randi([0 1], 1, 64);
%% Initialization and Time Vector Generation
stream = [ 0 1];  % Use only 5 bits for clean waveform
bit_rate=1e3;
samples_per_rate=1000;  % Each bit will be represented by 100 samples in the time domain

num_bit=length(stream);
total_samples=num_bit * samples_per_rate;

t = linspace(0, num_bit / bit_rate, total_samples);

%% polar_nrz line coding  
polar_nrz=zeros(1,total_samples);

for i =1:num_bit
    start_index=(i-1)*samples_per_rate +1;
    end_index=i*samples_per_rate;
    polar_nrz(start_index:end_index) = 2 * stream(i) - 1;  % 1 -> +1, 0 -> -1
end

figure;
plot(t,polar_nrz, 'LineWidth', 1.5);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('polar NRZ Signal');
ylim([-1.2 1.2]);             
grid on;    

%%%%%%%%%%%%%%%%%%%%%%%BPSK%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implement a digital modulation system BPSK transmitter

%Use bitstream as your baseband input.:(polar_nrz line coding)
energy_bit= 1;
Tb= 1/bit_rate;
%A=sqrt((2*energy_bit)/Tb);
A=1;
fc = 10*bit_rate;
tb = (0:samples_per_rate-1) * (Tb / samples_per_rate);

signal_tx = zeros(1, total_samples);  
carrier_wave = A * cos(2 * pi * fc * tb);
for i = 1:num_bit
    start_index=(i-1)*samples_per_rate +1;
    end_index=i*samples_per_rate;
    if polar_nrz(start_index) == 1
        signal_tx(start_index:end_index) = carrier_wave;
    else
        signal_tx(start_index:end_index) = -carrier_wave;
    end
end
figure;
plot(t,signal_tx, 'LineWidth', 1);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('BPSK Transmitted Signal');
ylim([-1.2*A 1.2*A]);                          
grid on;

%% Spectrum of BPSK Transmitted Signal
T=num_bit /bit_rate ;
fs = samples_per_rate*bit_rate;
ts= 1/fs;
df = fs/total_samples;
N=ceil (T/ts);

if mod(N, 2) == 0
    f=-(fs/2) : df : (fs/2)- df ;
else
    f=-((fs/2)-(df/2)) : df : (fs/2)- (df/2) ;
end

figure;
signal_fft = fftshift(fft(signal_tx))*(1/N);
plot(f/1000, abs(signal_fft));
xlabel('Frequency (kHz)');
ylabel('Magnitude');
title('Spectrum of BPSK Transmitted Signal');
grid on;
xlim([-20 20]);

%% Digital modulation system BPSK receiver 

% Create matching carrier wave for full duration
carrier_wave_rx = A * cos(2 * pi * fc * t); % Full-length carrier
v = signal_tx .* carrier_wave_rx;           % Coherent demodulation

% Initialize
rx_bits = zeros(1, num_bit);
decision_points = zeros(1, num_bit);

for i = 1:num_bit
    start_index = (i-1)*samples_per_rate + 1;
    end_index = i*samples_per_rate;

    v_segment = v(start_index:end_index);

    t_segment = t(start_index:end_index);  % Time vector for current bit
    integration_result = trapz(t_segment, v_segment);

    decision_points(i) = integration_result;

    % Decision based on sign
    if integration_result > 0
        rx_bits(i) = 1;
    else
        rx_bits(i) = -1;
    end
end

% Compare transmitted vs received bits
disp('Transmitted bits:');
disp(stream);
disp('Received bits:');
disp(rx_bits);

% Plot received signal decision points
figure;
stem(1:num_bit, decision_points, 'filled');
xlabel('Bit index');
ylabel('Integration result');
title('Receiver Decision Metric per Bit');
grid on;

rx_wave = zeros(1, total_samples);
for i = 1:num_bit
    start_index = (i-1)*samples_per_rate + 1;
    end_index = i*samples_per_rate;
    rx_wave(start_index:end_index) = rx_bits(i);
end

figure;
plot(t, rx_wave, 'LineWidth', 1.5);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('BPSK Received Signal (Waveform)');
ylim([-1.2 1.2]);
grid on;

%% Spectrum of BPSK Transmitted Signal
T=num_bit /bit_rate ;
fs = samples_per_rate*bit_rate;
ts= 1/fs;
df = fs/total_samples;
N=ceil (T/ts);

if mod(N, 2) == 0
    f=-(fs/2) : df : (fs/2)- df ;
else
    f=-((fs/2)-(df/2)) : df : (fs/2)- (df/2) ;
end

figure;
rx_wave_fft = fftshift(fft(rx_wave))*(1/N);
plot(f/1000, abs(rx_wave_fft));
xlabel('Frequency (kHz)');
ylabel('Magnitude');
title('Spectrum of BPSK Received Signal');
grid on;
xlim([-20 20]);

%% Digital modulation system BPSK receiver with phases
phase = [30 60 90] * pi/180;

% Initialize
rx_bits = zeros(1, num_bit);
decision_points = zeros(1, num_bit);
for j = 1:3
    ph = phase(j);
    carrier_wave_rx = A * cos((2 * pi * fc * t)+ph); % Full-length carrier
v = signal_tx .* carrier_wave_rx;           % Coherent demodulation
for i = 1:num_bit
    start_index = (i-1)*samples_per_rate + 1;
    end_index = i*samples_per_rate;

    v_segment = v(start_index:end_index);

    t_segment = t(start_index:end_index);  % Time vector for current bit
    integration_result = trapz(t_segment, v_segment);

    decision_points(i) = integration_result;

    % Decision based on sign
    if integration_result > 0
        rx_bits(i) = 1;
    else
        rx_bits(i) = -1;
    end
end
    % Compare transmitted vs received bits
    disp('Transmitted bits:');
    disp(stream);
    fprintf('Received bits at phase %.0f° \n', ph * 180 / pi);
    disp(rx_bits);
    
    % Plot received signal decision points
    figure;
    stem(1:num_bit, decision_points, 'filled');
    xlabel('Bit index');
    ylabel('Integration result');
    title(sprintf('Receiver Decision Metric per Bit at phase %.0f°',  ph*180/pi));

    grid on;
    
    rx_wave = zeros(1, total_samples);
    for i = 1:num_bit
        start_index = (i-1)*samples_per_rate + 1;
        end_index = i*samples_per_rate;
        rx_wave(start_index:end_index) = rx_bits(i);
    end
    
    figure;
    plot(t, rx_wave, 'LineWidth', 1.5);
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    title(sprintf('BPSK Received Signal (Waveform) at phase %.0f°',  ph*180/pi));
    ylim([-1.2 1.2]);
    grid on;
    
    
    T=num_bit /bit_rate ;
    fs = samples_per_rate*bit_rate;
    ts= 1/fs;
    df = fs/total_samples;
    N=ceil (T/ts);
    
    if mod(N, 2) == 0
        f=-(fs/2) : df : (fs/2)- df ;
    else
        f=-((fs/2)-(df/2)) : df : (fs/2)- (df/2) ;
    end
    
    figure;
    rx_wave_fft = fftshift(fft(rx_wave))*(1/N);
    plot(f/1000, abs(rx_wave_fft));
    xlabel('Frequency (kHz)');
    ylabel('Magnitude');
    title('Spectrum of BPSK Received Signal');
    title(sprintf('Spectrum of BPSK Received Signal at phase %.0f°',  ph*180/pi));
    grid on;
    xlim([-20 20]);

    % Convert transmitted stream from [1 0 1 ...] to [+1 -1 +1 ...] for comparison
tx_bits = 2*stream - 1;

% Calculate number of bit errors
num_errors = sum(rx_bits ~= tx_bits);

% Bit Error Rate (BER)
ber = num_errors / num_bit;

fprintf('Number of bit errors at phase %.0f°: %d\n', ph * 180 / pi, num_errors);
fprintf('Bit Error Rate (BER) at phase %.0f°: %.4f\n', ph * 180 / pi, ber);

end
