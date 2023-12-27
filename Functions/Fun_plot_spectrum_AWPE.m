function Fun_plot_spectrum_AWPE(Sig_1, Sig_2, Sig_3, Fs)

x_max = ceil(length(Sig_3)/Fs); 
f_max = min(Fs/2/1e3,8);

N1 = length(Sig_1); % length of sig_1, in sampling points
N2 = length(Sig_2); % length of sig_2, in sampling points
N3 = length(Sig_3); % length of sig_2, in sampling points
time_1 = (0:N1-1)/Fs; % time of sig_1, in second
time_2 = (0:N2-1)/Fs; % time of sig_2, in second
time_3 = (0:N3-1)/Fs; % time of sig_2, in second

Barmin = -40; % minimum of color bar, in dB
Barmax = 30; % maximum of color bar, in dB

%% 3. generate the waveform and spectrum

r    = 10; % resoultion, in Hz
nfft = fix(Fs/r); % fft
win  = hamming(nfft); % window
olap = fix(nfft*0.5); % overlap

% generate the spectrum
[Y_1, freq, frameTime] = spectrogram(Sig_1, win, olap, 256, Fs); 
[Y_2, ~, ~]   = spectrogram(Sig_2, win, olap, 256, Fs);
[Y_3, ~, ~]   = spectrogram(Sig_3, win, olap, 256, Fs);

% 4-1. waveform of clean signal

linewd = 1.5;
figure('position',[50, 30, 1140, 550]);
subplot(321)
t_vect   = [1:length(Sig_1)]'/Fs;
plot(t_vect, Sig_1,'-b','LineWidth',linewd);grid on;
text(20,-0.7,'(a)','fontName','Times New Roman','FontSize',14,'HorizontalAlignment','center')
axis([0 x_max -0.5 0.5]);


subplot(322)
imagesc(frameTime,freq./1e3,20*log10(abs(Y_1))); %
axis xy; 
ylabel('$f~(\mathrm{kHz})$', 'interpreter', 'latex');
axis([0 x_max 0 f_max]);
text(20,-1,'(b)','fontName','Times New Roman','FontSize',14,'HorizontalAlignment','center')
set(gca, 'fontSize', 12);
set(gca, 'fontName', 'Times New Roman');
colormap jet;
caxis([Barmin Barmax]);

subplot(323)
plot(t_vect, Sig_2,'-b','LineWidth',linewd);grid on;
text(20,-0.7,'(c)','fontName','Times New Roman','FontSize',14,'HorizontalAlignment','center')
axis([0 x_max -0.5 0.5]);

subplot(324)
imagesc(frameTime,freq./1e3,20*log10(abs(Y_2))); % 
axis xy; 
ylabel('$f~(\mathrm{kHz})$', 'interpreter', 'latex');
xlabel('Time~$(\mathrm{s})$', 'interpreter', 'latex','position',[5,-1.5,0]);
axis([0 x_max 0 f_max]);
text(20,-1,'(d)','fontName','Times New Roman','FontSize',14,'HorizontalAlignment','center')
set(gca, 'fontSize', 12);
set(gca, 'fontName', 'Times New Roman');
caxis([Barmin Barmax]);

subplot(325)
plot(t_vect, Sig_3,'-b','LineWidth',linewd);grid on;
text(20,-0.7,'(e)','fontName','Times New Roman','FontSize',14,'HorizontalAlignment','center')
axis([0 x_max -0.5 0.5]);

subplot(326)
imagesc(frameTime,freq./1e3,20*log10(abs(Y_3))); % 
axis xy; 
ylabel('$f~(\mathrm{kHz})$', 'interpreter', 'latex');
xlabel('Time~$(\mathrm{s})$', 'interpreter', 'latex','position',[5,-1.5,0]);
axis([0 x_max 0 f_max]);
text(20,-1,'(f)','fontName','Times New Roman','FontSize',14,'HorizontalAlignment','center')
set(gca, 'fontSize', 12);
set(gca, 'fontName', 'Times New Roman');
% colormap jet;colorbar('south','position',[0.09 0.08 0.83 0.05]);
caxis([Barmin Barmax]);



