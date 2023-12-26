function [Xsignal, Xclean, Xnoise] = add_Multinoise(clean_signal, noise_signal, SNR)
% clean_signal£º  clean speech signal
% noise_signal:   noise signal
% SNR: Signal Noise Ratio

% Authors: Gongping Huang
% Date: 11/11/2013
%-------------------------------------------------------------------------
Xclean = clean_signal;
noise = noise_signal;
xLen = length(Xclean);

sigmax2 = Xclean(:,1)'*Xclean(:,1)/xLen;
sigmav2 = noise(:,1)'*noise(:,1)/xLen;
scale = sqrt(sigmax2 * 10^(-SNR/10)/sigmav2);
Xnoise = scale* noise;  % The noise power at each channel should always the same
Xsignal = Xclean + Xnoise;
end



