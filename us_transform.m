%% Processing ultrasound signal with fft and wavelet transform
% preliminary test for pilot study for DFG grant
% input: ultrasound A-mode signal

clear all;
clear figure;

filename='sensor14.xlsx';
A=xlsread(filename);
N=size(A,1);
Data=A(:,2);
US=A(:,2);
time=A(:,1);
plot(A(:,1),A(:,2),'DisplayName','A(:,1:2)'), hold on;

%% Signal smoothing - low pass filter butterworth
% 
freq = 2e06; 
% [b,a]=butter(2,70000/freq,'low');  
% [b,a]=butter(2,100000/freq,'high'); 
[b,a]=butter(2,[60000 1500000]/freq); 
FD=filtfilt(b,a,US);         % filtered data

figure(2);
plot(US,'y'); hold on; plot(FD); hold off; ylim([-1,1]); title('Low pass filtered');

%% Hilbert Transform
% envelope(signal,Fs);
% HD=hilbert(FD);
HD=envelope(FD,freq);
HilHD=abs(HD);          

figure(3);
plot(time,HilHD); ylim([0 1]); xlabel('time/s'); xlim([0 1.6e-04]); ylabel('amplitude/V'); title('Hilbert-transformed data');

%% STFT

D_f=HilHD;
window = 50;        % epoch is 0.5 sec. Hamming window of length 1024.   

% spectrogram(D_f,window,window/2,window,freq,'yaxis'); ylim([0 1e5]); 
[ST,FT,TI,PS] = spectrogram(US,window,window/2,window,freq); %evaluate the spectrum at [1024/2+1]=513 freq
    % ST - short-time fourier transform, FT - cyclical freq wrt to fs, 
    % TI - time instants vector, PS - PSD estimate 
    
figure(4);
subplot(2,1,1); plot(time(:,1),US(:,1)); ylim([-1 1]); xlim([0 1.6e-04]); %ylim([0 1.0]); ylabel('Hilbert'); 
subplot(2,1,2); spectrogram(US,window,window/2,window,freq,'yaxis'); %ylim([0 1e6]); 
colormap bone

%% CWT
coefs = cwt(US,1:100,'sym2');
figure(5);
subplot(2,1,1); plot(time(:,1),HilHD(:,1)); ylabel('Hilbert'); ylim([0 1.0]); xlim([0 1.6e-04]);
subplot(2,1,2); cwt(US,1:100,'sym2','scal'); xlim([0 16000]); title('Scalogram'); ylabel('Scale')

%% Fast fourier transform

Y= fft(A);
t = linspace(0,1,N);
T = t(2)-t(1);
Fs = 1/T;
f = Fs/2*linspace(0,1,N);
figure(6);
plot(f, abs(Y));

%% 
% Fs=1000;
% T = 1/Fs;   
% t = (0:len-1)*T; 
NFFT = 2^nextpow2(len); % Next power of 2 from length of y
Y = fft(A1,NFFT)/len;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

