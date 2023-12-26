function [z, delay] = derev_awpe_RLS(y,Para_AWPE)

% implementation the recursive least square (RLS) based adaptive weighted
% prediction error (AWPE) dereverberation algorithm for speech enhancement.
%
% In this implementation, we divide the speech enhancement system into
% four steps:
% 1) Pre-beaforming: we use a fixed beamforming--superdirective beamforming
% 2) Pre-noise reduction: we use parameteric Wiener filter
% 3) Dereverberation: we use recursive least square (RLS) based AWPE
% 4) Post-noise reduction: we use another parameteric Wiener filter
%
% USAGE:
% [z, delay] = derev_wpe_RLS(y,Para_method,Para_AWPE,Para_FBF,Para_NR)
%
% Inputs summary:
% y:            the mulichannel noisy signals
% Para_method:  parameters for method choosing
% Para_AWPE:    parameters for AWPE
% Para_FBF:     parameters for fixed beamforming
% Para_NR:      parameters for pre-/post-noise reduction
%
% some important parameters in AWPE
% Kl:           the order of the prediction filter
% delta_N:      the delay in frame
% alpha_rls:    the forgetting factor

% ----------------------------------------------------------------
% Output:
% z        the enhanced signal
% delay:   time delay
% ----------------------------------------------------------------
% Notice:
%
% 1) By introducing a certain delay,namely,delta_N to prevent the
% STFTs of the current and past observed signals from sharing
% the same signal in the time domain. Because of this delay,
% the algorithm cannot predict early reflections of reverberation.
%
% 2) We use pre-beaforming, pre-noise reduction, and post-noise reduction to
% deal with the noisy signal, it may be adjusted depend on the actual
% environment, i.e, if we need a noise reduction. The noise PSD is
% estimated from the noisy observation and it used to construct the
% Wiener filter, where the two Wiener filter are used in the beamforming
% output and dereverberated signal, respectively.
%
% 3) we may set the order of the prediction filter varied with frequency
% as the author suggested in [1], i.e., a higher value at low frequencies
% and lower value at high frequencies.

% ----------------------------------------------------------------
% Reference,
% the AWPE part:
% [1] T Yoshioka, ��Speech enhancement in reverberant environments��, Ph.D.
% dissertation, Kyoto University, 2010
% the FBF part:
% [2] G. Huang, J. Benesty, and J. Chen, "Subspace superdirective
% beamforming with uniform circular microphone arrays," in Proc.
% IEEE IWAENC, pp. 1-5, 2016.
% the NR part:
% [3] I. Cohen, "Noise spectrum estimation in adverse environments:
% 	Improved minima controlled recursive averaging," IEEE Transactions
% on Speech and Audio Processing, 11(5), 466-475, 2003.
%
% Authors: Gongping Huang
% Date: 07/30/2017
%-------------------------------------------------------------------------
%

%% % 1) Get some parameters

% 1-1)parameters for method choosing
% do_iteration        = Para_AWPE.do_iteration; % if do iteration
% 1-2)some fundamental parameters
M                   = length(y(1,:));
% 1-3)parameters for STFT
% frame and window sizes
frmsize             = 64*1;
noverlap            = 4;
winsize             = frmsize*noverlap;
nfbands             = winsize/2+1;
delay               = (noverlap-1)*frmsize;
% Window parameters
kwin                = kaiser(winsize, 1.9*pi);  % Kaiser Window
kwsigma             = sqrt(sum(kwin.^2)/winsize*noverlap);
kwin                = kwin/kwsigma;

% 1-5)parameters for AWPE
Kl1                 = Para_AWPE.Kl;
delta_N             = Para_AWPE.delta_N;
alpha_rls           = Para_AWPE.alpha_rls;
N                   = Kl1 + delta_N;
lamuda_min          = Para_AWPE.lamuda_min;

inv_alpha_rls = 1/alpha_rls;
% 1-7)reconstruct the noisy signal-just for simulation
K                   = length(y(:,1));
nfrms               = ceil(K/frmsize);
Kfrm                = nfrms*frmsize;
if (Kfrm > K)
    y(K+1:Kfrm,:)   = 0;
end

% 1-8)Create the needed buffers for FFT/IFFT
inBuf               = zeros(winsize,M);   % input buffer
outWin              = zeros(winsize,1);   % output analysis window
outOAWin            = zeros(winsize+frmsize,1); % output overlap-add window
inWinBuf            = zeros(M*N,winsize);   % input fft buffer

% set some initial some parameters for dereverberation
% initalize the filter vector with zero vector
h_matrix1           = zeros(M*Kl1,nfbands);
% initalize the convariance matrix with identity matrix
Cov_y1              = zeros(M*Kl1,M*Kl1,nfbands);
for nbin= 1 :nfbands
    Cov_y1(:,:,nbin)  = eye(M*Kl1);
end

% 1-10)Create the needed buffers for finally output
z                   = zeros(Kfrm,1);

% to show the process time
waitHandle = waitbar(0,'Please wait...');

% re-initalize the convariance matrix with identity matrix
Cov_y1              = zeros(M*Kl1,M*Kl1,nfbands);
for nbin= 1 :nfbands
    Cov_y1(:,:,nbin)  = eye(M*Kl1);
end

% %---------------------------------------------------------------
%% % speech enhancement process
tic;
for nf = 1:nfrms
    for m =  M : -1 : 1
        yfrm                = y((nf-1)*frmsize+1:nf*frmsize,m); % get a frame of data
        inBuf(:,m)          = [inBuf(frmsize+1:end,m); yfrm];   % shift it in
        inWin               = fft(inBuf(:,m).*kwin);        % apply the kaiser window
        inWinBuf(2:end,:)   = inWinBuf(1:end-1,:);
        inWinBuf(1,:)       = inWin;
    end
    
    for nbin = 1 :nfbands %winsize        
        Fyvect          = inWinBuf(M*delta_N+1:end,nbin);
        Xr              = inWinBuf(1,nbin);
        
        % estimate the speech model parameters
        s_out           = Xr - h_matrix1(:,nbin)'*Fyvect;
        
        % calculate the PSD, here we use a smotthing and exponential factor
        % to improve the performance
        lamuda       = abs(s_out)^2;
        %lamuda       = max(lamuda,lamuda_min);
        
        % update the weighted cross-/correlation matrix
        Bn           = Cov_y1(:,:,nbin)*Fyvect;
        B            = 1/max(alpha_rls*lamuda + Fyvect'*Bn,1e-6);
        Kf           = B.*Bn;
        % update the covariance matrix
        An           = Fyvect'*Cov_y1(:,:,nbin);
        % An           = Bn';
        Cov_y1(:,:,nbin)  = inv_alpha_rls*(Cov_y1(:,:,nbin) - Kf*An);
        % update the mean vector
        h_matrix1(:,nbin) = h_matrix1(:,nbin) + Kf*conj(s_out);
        
        % re-implement the dereverberation (can improve the performance)
        % s_out           = Xr - h_matrix1(:,nbin)'*Fyvect;
        outWin(nbin)      = s_out;
        
    end
    outWin(nfbands+1:end)   = conj(flipud(outWin(2:nfbands-1)));
    % inverse fft
    outWin = real(ifft(outWin));
    % overlap add
    outOAWin(frmsize+1:end) = outOAWin(frmsize+1:end)+ kwin.*outWin;
    outOAWin(1:winsize)     = outOAWin(frmsize+1:end);
    outOAWin(winsize+1:end) = 0;
    % copy to the output buffer
    outBuf  = (outOAWin(1:frmsize));
    z((nf-1)*frmsize+1 : nf*frmsize) = outBuf;
    waitbar(nf/nfrms);
end
z = z(:);
close(waitHandle)
toc
end
