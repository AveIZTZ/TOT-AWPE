function [h_matrix1, Cov_y1, h_matrix2, Cov_y2] = DR_initalize_kron(y, Para_STFT, Para_AWPE)

% Initalize the AWPE filters

M           = length(y(1,:));
% parameters for STFT
frmsize     = Para_STFT.frmsize;
winsize     = Para_STFT.winsize;
nfbands     = Para_STFT.nfbands;
kwin        = Para_STFT.kwin;

% parameters for AWPE
Fs           = Para_AWPE.Fs;
Kl1          = Para_AWPE.Kl;
delta_N      = Para_AWPE.delta_N;
alpha_rls    = Para_AWPE.alpha_rls;
N            = Kl1 + delta_N;


% Italization for Kronecker filters
alpha       = alpha_rls;
P           = Para_AWPE.P;
L1          = Para_AWPE.L1;
L2          = Para_AWPE.L2;

% initalize the dereverberation filters with zero vector
h_matrix1   = zeros(P*L1,nfbands);
h_matrix2   = zeros(P*L2,nfbands);
% initalize the convariance matrix with identity matrix
Cov_y1      = zeros(P*L1,P*L1,nfbands);
Cov_y2      = zeros(P*L2,P*L2,nfbands);

% with Kronecker filter, hhow to initalize the filters is important
% may need future optimization
ct          = Para_AWPE.ct;
h2p         = [ct*eye(P);zeros(L2-P,P)];
% h2p         = [ct*ones(P);zeros(L2-P,P)];

for nbin= 1 :nfbands
    Cov_y1(:,:,nbin)  = eye(P*L1);
    Cov_y2(:,:,nbin)  = eye(P*L2);
    h_matrix2(:,nbin) = h2p(:);
end

xk_1         = zeros(P*L1,1);
xk_2         = zeros(P*L2,1);
IM_1         = eye(L1);
IM_2         = eye(L2);

% reconstruct the noisy signal
K            = length(y(:,1));
if K < winsize
    fprintf('Without initialization! \n');
    return
else
    fprintf('With initialization! \n');
end

nfrms               = ceil(K/frmsize);
Kfrm                = nfrms*frmsize;
if (Kfrm > K)
    y(K+1:Kfrm,:)   = 0;
end

% Create the buffers for FFT/IFFT
inBuf          = zeros(winsize,M);   % input buffer
inWinBuf       = zeros(M*N,winsize);   % input fft buffer

% recursive factor
inv_alpha_rls  = 1/alpha_rls;
lamuda_min     = 1e-5;


for nframe = 1:nfrms
    
    
    for m =  M : -1 : 1
        yfrm                = y((nframe-1)*frmsize+1:nframe*frmsize,m); % get a frame of data
        inBuf(:,m)          = [inBuf(frmsize+1:end,m); yfrm];   % shift it in
        inWin               = fft(inBuf(:,m).*kwin);        % apply the kaiser window
        inWinBuf(2:end,:)   = inWinBuf(1:end-1,:);
        inWinBuf(1,:)       = inWin;
    end
    
    
    for nbin = 1 :nfbands %winsize
        % get the n1th row of the buffer
        Fyvect          = inWinBuf(M*delta_N+1:end,nbin);
        % construct the delayed vector
        Xr              = inWinBuf(1,nbin);
        
        % update as a 2-way Tensor/Kronecker
        % 1) update for n = 1
        % 1-1) compute Hk1_p and xk_1
%         for p=1:P
%             Hk1_p                   = kron(h_matrix2((p-1)*L2+1:p*L2,nbin),IM_1);
%             xk_1((p-1)*L1+1:p*L1)   = Hk1_p'*Fyvect;
%         end
        Y = reshape(Fyvect,[L1,L2]);
        for p=1:P
            xk_1((p-1)*L1+1:p*L1)   = Y*conj(h_matrix2((p-1)*L2+1:p*L2,nbin));
        end
        % Gusee we don't need to compute two correlation matrix,
        % this should be simpled! -- June 10, 2021
        
        % Compute the output signal
        s_out1        = Xr - h_matrix1(:,nbin)'*xk_1;
        conj_s_out1   = conj(s_out1);
                
        % calculate the PSD
        lamuda1       = max(s_out1*conj_s_out1,lamuda_min);
        % update the weighted cross-/correlation matrix
        Bn_1         = Cov_y1(:,:,nbin) * xk_1;
        An_1         = xk_1'*Cov_y1(:,:,nbin);
        An_d         = 1/max(alpha*lamuda1 + xk_1'*Bn_1, 1e-6);
        Kf_1         = An_d.*Bn_1;
        
        Cov_y1(:,:,nbin)     = inv_alpha_rls.*(Cov_y1(:,:,nbin) - Kf_1*An_1);
        % update the filter
        h_matrix1(:,nbin)        = h_matrix1(:,nbin) + conj_s_out1.*Kf_1;
        
        % 2) update for n = 2
        % 2-1) compute H2_p and xk_2
%         for p=1:P
%             Hk2_p                   = kron(IM_2,h_matrix1((p-1)*L1+1:p*L1,nbin));
%             xk_2((p-1)*L2+1:p*L2)   = Hk2_p'*Fyvect;
%         end
        for p=1:P
            xk_2p   = h_matrix1((p-1)*L1+1:p*L1,nbin)'*Y;
            xk_2((p-1)*L2+1:p*L2)   = xk_2p(:);
        end
        % 2-2) Compute the output signal
        s_out2       = Xr - h_matrix2(:,nbin)'*xk_2;
        conj_s_out2  = conj(s_out2);
        
        % calculate the PSD
        lamuda2      = max(s_out2*conj_s_out2,lamuda_min);
        
        % update the weighted cross-/correlation matrix
        Bn_2         = Cov_y2(:,:,nbin) * xk_2;
        An_2         = xk_2'*Cov_y2(:,:,nbin);
        Bn_d         = 1/max(alpha*lamuda2 + xk_2'*Bn_2, 1e-6);
        Kf_2         = Bn_d.*Bn_2;
        
        Cov_y2(:,:,nbin)     = inv_alpha_rls.*(Cov_y2(:,:,nbin) - Kf_2*An_2);
        % update the filter
        h_matrix2(:,nbin)    = h_matrix2(:,nbin) + conj_s_out2.*Kf_2;
        
    end % for nbin
end

end
