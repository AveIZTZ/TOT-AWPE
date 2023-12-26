function x = simosys(s, H, M, L)
% SIMOSYS: Single-Input Multiple-Output FIR System
%
% Usage: x = simosys(s, H, M, L)
%   Input:
%     s : the source signal.
%     H : the multi-channel impulse responses (in the time domain).
%           H = [h_1 h_2 ... h_M]
%     M : the number of channels.
%     L : the length of a channel impulse response.
%   Output:
%     x : the output of the SIMO FIR system.
%
% ------
% Author:  Jingdong Chen, Ph.D.
%          Member of Technical Staff
%          Bell Labs, Lucent Technologies
% Copyright 2002(c)
%
if (nargin ~= 4)
    error('Incorrect number of inputs!');
end
if (M <= 0)
    error('Incorrect number of channels!');
end
if (L <= 0)
    error('Incorrect length of a channel impulse response!');
end
[R, C] = size(H);
if (R ~= L)
    error('The row number of H doesn''t match L!');
end
if (C ~= M)
    error('The column number of H doesn''t match M!');
end

L2 = 2*L;

s  = s(:);
Ndata = length(s);
if (Ndata <= L2)
    for m = 1:M
        x(:,m) = filter(H(:,m),1,s);
    end
elseif L <= 4
    for m = 1:M,
        x(:,m) = filter(H(:,m),1,s);
    end;
else
    HF = fft([H; zeros(L,M)]);
    Nblocks = ceil(Ndata/L);
    sb = zeros(L2,1);
    for n = 1:Nblocks
        if (n ~= Nblocks)
            sb = [sb(L+1:end); s((n-1)*L+1:n*L)];
        else
            ndleft = Ndata - (Nblocks-1)*L;
            sb = [sb(L+1:end); zeros(L,1)];
            sb(L+1:L+ndleft) = s((n-1)*L+1:end);
        end
        SBF = fft(sb);
        for m = 1:M
            XBF(:,m) = SBF.*HF(:,m);  % 使用的是重叠保留法
        end
        xb  = real(ifft(XBF));
        if (n == 1)
            x = xb(L+1:end,:);
        else
            x = [x; xb(L+1:end,:)];
        end
    end
end
x = x(1:Ndata, :);  % 建立过程考虑，因为数据无限长，所以不考虑结束过程
        