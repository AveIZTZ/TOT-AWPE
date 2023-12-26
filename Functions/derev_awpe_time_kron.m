function [z, delay] = derev_awpe_time_kron(y, Para_AWPE)
% implementation the recursive least square (RLS) based Kronecker adaptive weighted
% prediction error (AWPE) dereverberation algorithm for speech enhancement.
%
% ----------------------------------------------------------------
%
% Authors: Gongping Huang
% Date: 02/03/2022
%-------------------------------------------------------------------------
%

% 1-1) load the speaker change information
M                  = length(y(1,:));

% 1-3) parameters for STFT
frmsize             = 128*1; % frame size
noverlap            = 4;
winsize             = frmsize*noverlap; % window size
nfbands             = winsize/2+1;
delay               = (noverlap-1)*frmsize;
% Window parameters
kwin                = kaiser(winsize, 1.9*pi);  % Kaiser Window
kwsigma             = sqrt(sum(kwin.^2)/winsize*noverlap);
kwin                = kwin/kwsigma;

% 1-5) parameters for AWPE
Fs                  = Para_AWPE.Fs;
Kl1                 = Para_AWPE.Kl;
delta_N             = Para_AWPE.delta_N;
alpha_rls           = Para_AWPE.alpha_rls;
N                   = Kl1 + delta_N;
% minimum value of the variance on the denominator
lamuda_min          = Para_AWPE.lamuda_min;

% 1-6) reconstruct the noisy signal-just for simulation
K                   = length(y(:,1));
nfrms               = ceil(K/frmsize);
Kfrm                = nfrms*frmsize;
if (Kfrm > K)
    y(K+1:Kfrm,:)   = 0;
end

% 1-7) Create the needed buffers for finally output
z                   = zeros(Kfrm,1);

% 1-9) inverse of the recursive factor
inv_alpha_rls_start = 1/alpha_rls;
inv_alpha_rls       = inv_alpha_rls_start;

% Italization for Kronecker filters
alpha               = alpha_rls;
P                   = Para_AWPE.P;
L1                  = Para_AWPE.L1;
L2                  = Para_AWPE.L2;

xk_1                = zeros(P*L1,1);
xk_2                = zeros(P*L2,1);

% % ##################################################################
% % # %%            Start dereverberation process                    #
% % ##################################################################

% to show the process time
waitHandle          = waitbar(0,'Please wait...');

inBuf               = zeros(winsize,M);         % input buffer
outWin              = zeros(winsize,1);         % output analysis window
outOAWin            = zeros(winsize+frmsize,1); % output overlap-add window
inWinBuf            = zeros(M*N,winsize);       % input fft buffer

% initalize the dereverberation filters with zero vector
h_matrix_1   = zeros(P*L1,nfbands);
h_matrix_2   = zeros(P*L2,nfbands);
% initalize the convariance matrix with identity matrix
Cov_y1      = zeros(P*L1,P*L1,nfbands);
Cov_y2      = zeros(P*L2,P*L2,nfbands);
% with Kronecker filter, hhow to initalize the filters is important, may need future optimization
ct          = Para_AWPE.ct;
h2p         = [ct*eye(P);zeros(L2-P,P)];
for nbin= 1 :nfbands
    Cov_y1(:,:,nbin)  = eye(P*L1);
    Cov_y2(:,:,nbin)  = eye(P*L2);
    h_matrix_2(:,nbin) = h2p(:);
end


tic;

for nt = 1:length(y1)
    
    
    for m =  M : -1 : 1
        yfrm                = y((nf-1)*frmsize+1:nf*frmsize,m); % get a frame of data
        inBuf(:,m)          = [inBuf(frmsize+1:end,m); yfrm];   % shift it in
        inWin               = fft(inBuf(:,m).*kwin);            % apply the kaiser window
        inWinBuf(2:end,:)   = inWinBuf(1:end-1,:);
        inWinBuf(1,:)       = inWin;
    end
        
    
    for nbin = 1 :nfbands %winsize
        % get the n1th row of the buffer
        Fyvect          = inWinBuf(M*delta_N+1:end,nbin);
        % construct the delayed vector
        Xr              = inWinBuf(1,nbin);
        
        % 1) update g1
        Y = reshape(Fyvect,[L1,L2]);
        for p=1:P
            xk_1((p-1)*L1+1:p*L1)   = Y*conj(h_matrix_2((p-1)*L2+1:p*L2,nbin));
        end
        
        % Compute the output signal
        s_out1        = Xr - h_matrix_1(:,nbin)'*xk_1;
        conj_s_out1   = conj(s_out1);
        
        % calculate the PSD
        % lamuda1       = max(s_out1*conj_s_out1,lamuda_min);
        lamuda1       = s_out1*conj_s_out1;
        
        % update the weighted cross-/correlation matrix
        Bn_1          = Cov_y1(:,:,nbin) * xk_1;
        An_1          = xk_1'*Cov_y1(:,:,nbin);
        An_d          = 1/max(alpha*lamuda1 + xk_1'*Bn_1, 1e-6);
        Kf_1          = An_d.*Bn_1;
        
        Cov_y1(:,:,nbin)     = inv_alpha_rls.*(Cov_y1(:,:,nbin) - Kf_1*An_1);
        % update the filter
        h_matrix_1(:,nbin)    = h_matrix_1(:,nbin) + conj_s_out1.*Kf_1;
        
        
        % 2) update g2
        for p=1:P
            xk_2p   = h_matrix_1((p-1)*L1+1:p*L1,nbin)'*Y;
            xk_2((p-1)*L2+1:p*L2)   = xk_2p(:);
        end
        
        % 2-2) Compute the output signal
        s_out2       = Xr - h_matrix_2(:,nbin)'*xk_2;
        conj_s_out2  = conj(s_out2);
        
        % calculate the PSD
        % lamuda2      = max(s_out2*conj_s_out2,lamuda_min);
        lamuda2      = s_out2*conj_s_out2;

        % update the weighted cross-/correlation matrix
        Bn_2         = Cov_y2(:,:,nbin) * xk_2;
        An_2         = xk_2'*Cov_y2(:,:,nbin);
        Bn_d         = 1/max(alpha*lamuda2 + xk_2'*Bn_2,1e-6);
        Kf_2         = Bn_d.*Bn_2;
        
        Cov_y2(:,:,nbin)     = inv_alpha_rls.*(Cov_y2(:,:,nbin) - Kf_2*An_2);
        % update the filter
        h_matrix_2(:,nbin)    = h_matrix_2(:,nbin) + conj_s_out2.*Kf_2;
        
        % re-implement the dereverberation can improve the performance
        outWin(nbin)       = s_out2;
        
    end % for nbin
    
    % % ##################################################################
    % % # %%            Synthesis with overlap save                      #
    % % ##################################################################
    
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