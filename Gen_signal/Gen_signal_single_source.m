clear; clc;
close all;

load('impulse_response_0_180.mat');

L_early  = 0.04*Fs;
H1_early = H1(1:L_early,:);% source location -- 0 degree
figure('position',[50, 30, 1140, 550]);
subplot(211);plot(H1);
subplot(212);plot(H1_early);

saveXsignal         = 1;
% find the delay of direct path 
[~,D_delay1] = max(H1(:,1)); % D_delay1 = 92;

% % load the clean signal
Lentime1     = 60;       % the wave length which will be filtered
PathWave1    = 'women_16k.wav';% the sound wave which is the clean signal
[x1, fs0]    = audioread(PathWave1);% read the wave file
LS           = Lentime1*fs0;% calculate the signal length in samples
Sclean      = resample(x1(1 :LS), Fs, fs0);      %resample the clean signal with sample rate fs
 
% Generate the clean signal and reference signal with early reflections (40 ms) as 5 segments
[rowH1, cloumnH1] = size(H1);
[rowHE1, cloumnHE1] = size(H1_early);


Xclean = simosys(Sclean, H1, cloumnH1, rowH1);
Xearly = simosys(Sclean, H1_early, cloumnHE1, rowHE1);
% align the signal
Xearly      = Xearly(D_delay1:end,:); 
Xclean      = Xclean(D_delay1:end,:); 
Sclean      = Sclean(1:length(Xclean(:,1)));
 
% % show the generated reverberant signal
t = [1:length(Sclean(:,1))]/Fs;
figure('position',[0 50 600 400]);linewd=1.25;
plot(t,Xclean(:,1),'-r','LineWidth',linewd);hold on;grid on;
plot(t,Sclean(:,1),'-b','LineWidth',linewd);

if saveXsignal
    if Fs == 8000
    save('cleansignal_8k.mat', 'Fs','Sclean','Xclean','Xearly','rcoeff','delta');
    else
    save('cleansignal_16k.mat', 'Fs','Sclean','Xclean','Xearly','rcoeff','delta');
    end
end
