function mx=gcc(x,y,gama,L,ro,maxtd)
% GCC: the GCC Algorithm for Time Delay Estimation.
%
% Usage: mx=phat(x,y,gama,L,ro)
%
% mx: vector of delay estimate.
% x, y: two observation signals.
% gama: forgeting factor.
% L: frame length;
% ro: 0 for CC, and 1 for PHAT
% maxtd: maximum time delay. Delay estimate \in [-maxtd, maxtd]

if (length(x) ~= length(y))
  error('GCC: two signals should have the same length.');
end

if (ro ~=0 & ro ~= 1)
  error('GCC: ro has to be 0 for CC or 1 for PHAT.');
end

if (nargin == 5)
    maxtd = L/2-1;
end

SV = ones(L,1);
SV(maxtd+1:L-maxtd-1) = zeros(L-2*maxtd-1,1);    
n  = fix(size(x,1)/L);
Cxy = zeros(L,1);
nfbands  = L/2+1;
Hfilter = ones(L,1); Hfilter(1:L/16) = 0; %High-pass filter

for i = 1:n 
   % delay estimation 
   x1 = x(L*(i-1)+1:L*i).*kaiser(L,8);
   y1 = y(L*(i-1)+1:L*i).*kaiser(L,8);
   Xf = fft(x1);  Yf = fft(y1);
   Xf = Xf.*Hfilter;  Yf = Yf.*Hfilter;  %High-pass filtering
   Xf(nfbands+1:end) = conj(flipud(Xf(2:nfbands-1)));
   Yf(nfbands+1:end) = conj(flipud(Yf(2:nfbands-1)));
   Cxy = gama*Cxy + (1-gama)*(conj(Yf).*Xf);
   Eps = 1e-8;
   GCC = Cxy./(abs(Cxy).^ro + Eps);
   RR = real(ifft(GCC));
   RR = RR.*SV;
   [mv,mxx] = max(RR);
   if mxx < L/2,
      mx(i) = mxx - 1;
   else
      mx(i) = mxx - L-1;
   end
   if(mod(i,100)==0)
   if ro == 1,
       fprintf('\n PHAT: processing %d/%d',i,n);
   else
       fprintf('\n CC: processing %d/%d',i,n);
   end;
   end
end
