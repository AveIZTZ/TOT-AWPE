function [z, delay] = derev_awpe_RLS_TOT(y, Para_AWPE)
% implementation the recursive least square (RLS) based Kronecker adaptive weighted
% prediction error (AWPE) dereverberation algorithm for speech enhancement.
%
% ----------------------------------------------------------------
%
% Authors: Yujie Zhu
% Date: 11/08/2023
%-------------------------------------------------------------------------
%

% 1-1) load the speaker change information
M                  = length(y(1,:));

% 1-3) parameters for STFT
frmsize             = 64*1; % frame size
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
%lamuda_min          = Para_AWPE.lamuda_min;

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
L11                  = Para_AWPE.L11;
L12                  = Para_AWPE.L12;
L2                  = Para_AWPE.L2;

%xk_2_12                = zeros(P*L11*L2,1);
%xk_2_11                = zeros(P*L12*L2,1);
%xk_12_11                = zeros(L2*L2,1);

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
h_matrix_11   = zeros(L11*P*L2,nfbands);
h_matrix_12   = zeros(L12*P*L2,nfbands);
h_matrix_2   = zeros(L2*L2,nfbands);
% initalize the convariance matrix with identity matrix
Cov_y11      = zeros(P*L11*L2,P*L11*L2,nfbands);
Cov_y12      = zeros(P*L12*L2,P*L12*L2,nfbands);
Cov_y2      = zeros(L2*L2,L2*L2,nfbands);
% initalize the kronec product
%G2_12 = zeros(L11*L12*L2,P*L11*L2);
%G2_11 = zeros(L11*L12*L2,P*L12*L2);
%G12_11 = zeros(L11*L12*L2,L2*L2);
% with Kronecker filter, hhow to initalize the filters is important, may need future optimization
ct          = Para_AWPE.ct;

h12p = zeros(L12,P*L2);
h11p = zeros(L11,P*L2);
h2p = [ct*eye(L2)];
if L2 < L12
    for n = 1:L2
        h12p(n,(n-1)*P+1:n*P)         = ones(1, P);
    end
else
    for n = 1:L12
        h12p(n,(n-1)*P+1:n*P)         = ones(1, P);
    end
end
for nbin= 1 :nfbands
    Cov_y11(:,:,nbin)  = ct*eye(P*L11*L2);
    Cov_y12(:,:,nbin)  = ct*eye(P*L12*L2);
    Cov_y2(:,:,nbin)  = ct*eye(L2*L2);
    h_matrix_12(:,nbin) = h12p(:);
    h_matrix_11(:,nbin) = h11p(:);
    h_matrix_2(:,nbin) = h2p(:);
end
tic;

for nf = 1:nfrms
    
    % % ##################################################################
    % % # %%                      Implement STFT                         #
    % % ##################################################################
    
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

        % 1) update g2
        xk_12_11                = zeros(L2*L2,1);
        y_trans2              = reshape(Fyvect,[L11*L12,L2]);
        
        for k=1:L2
            G12_11_l = zeros(L11*L12,1);
            for q=1:P
                G12_11_lp = kron(h_matrix_12(((k-1)*P+q-1)*L12+1:((k-1)*P+q)*L12,nbin),h_matrix_11(((k-1)*P+q-1)*L11+1:((k-1)*P+q)*L11,nbin));
                G12_11_l = G12_11_l+G12_11_lp;
            end
            xk_12_11((k-1)*L2+1:k*L2) = conj(y_trans2'*G12_11_l);
        end
        
        
        % Compute the output signal
        s_out2        = Xr - h_matrix_2(:,nbin)'*xk_12_11;
        conj_s_out2   = conj(s_out2);
        
        % calculate the PSD
        %lamuda2       = max(s_out2*conj_s_out2,lamuda_min);
        lamuda2       = s_out2*conj_s_out2;
        
        % update the weighted cross-/correlation matrixh_matrix_2
        Bn_2          = Cov_y2(:,:,nbin) * xk_12_11;
        An_2          = xk_12_11'*Cov_y2(:,:,nbin);
        An_d          = 1/max(alpha*lamuda2 + xk_12_11'*Bn_2, 1e-6);
        Kf_2          = An_d.*Bn_2;
        
        Cov_y2(:,:,nbin)     = inv_alpha_rls.*(Cov_y2(:,:,nbin) - Kf_2*An_2);
        % update the filter
        h_matrix_2(:,nbin)    = h_matrix_2(:,nbin) + conj_s_out2.*Kf_2;      
                
        % 2) update g12
        xk_2_11                = zeros(P*L12*L2,1);
        y_trans12              = zeros(L11*L2,L12);
        for k=1:L2
            y_trans12((k-1)*L11+1:k*L11,:) = reshape(Fyvect((k-1)*L11*L12+1:k*L11*L12),[L11,L12]);
        end
        for k=1:L2
            for q=1:P
                G2_11_lp = kron(h_matrix_2((k-1)*L2+1:k*L2,nbin),h_matrix_11(((k-1)*P+q-1)*L11+1:((k-1)*P+q)*L11,nbin));
                xk_2_11(((k-1)*P+q-1)*L12+1:((k-1)*P+q)*L12) = conj(y_trans12'*G2_11_lp);
            end
        end
        
        % Compute the output signal
        s_out12        = Xr - h_matrix_12(:,nbin)'*xk_2_11;
        conj_s_out12   = conj(s_out12);
        
        % calculate the PSD
        %lamuda12       = max(s_out12*conj_s_out12,lamuda_min);
        lamuda12       = s_out12*conj_s_out12;
        
        % update the weighted cross-/correlation matrix
        Bn_12          = Cov_y12(:,:,nbin) * xk_2_11;
        An_12          = xk_2_11'*Cov_y12(:,:,nbin);
        An_d          = 1/max(alpha*lamuda12 + xk_2_11'*Bn_12, 1e-6);
        Kf_12          = An_d.*Bn_12;
        
        Cov_y12(:,:,nbin)     = inv_alpha_rls.*(Cov_y12(:,:,nbin) - Kf_12*An_12);
        % update the filter
        h_matrix_12(:,nbin)    = h_matrix_12(:,nbin) + conj_s_out12.*Kf_12;
       
        % 3) update g11
        xk_2_12                = zeros(P*L11*L2,1);
        y_trans11              = reshape(Fyvect,[L11,L12*L2]);
        for k=1:L2
            for q=1:P
                G2_12_lp = kron(h_matrix_2((k-1)*L2+1:k*L2,nbin),h_matrix_12(((k-1)*P+q-1)*L12+1:((k-1)*P+q)*L12,nbin));
                xk_2_12(((k-1)*P+q-1)*L11+1:((k-1)*P+q)*L11,:) = y_trans11*conj(G2_12_lp);
            end
        end
        
        % 2-2) Compute the output signal
        s_out11       = Xr - h_matrix_11(:,nbin)'*xk_2_12;
        conj_s_out11  = conj(s_out11);
        
        % calculate the PSD
        %lamuda11      = max(s_out11*conj_s_out11,lamuda_min);
        lamuda11      = s_out11*conj_s_out11;

        % update the weighted cross-/correlation matrix
        Bn_11         = Cov_y11(:,:,nbin) * xk_2_12;
        An_11         = xk_2_12'*Cov_y11(:,:,nbin);
        An_d         = 1/max(alpha*lamuda11 + xk_2_12'*Bn_11,1e-6);
        Kf_11         = An_d.*Bn_11;
        
        Cov_y11(:,:,nbin)     = inv_alpha_rls.*(Cov_y11(:,:,nbin) - Kf_11*An_11);
        % update the filter
        h_matrix_11(:,nbin)    = h_matrix_11(:,nbin) + conj_s_out11.*Kf_11;
        
        % re-implement the dereverberation can improve the performance
        outWin(nbin)       = s_out11;

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