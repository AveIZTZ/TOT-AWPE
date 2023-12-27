function [d, dn] = ItakuraSaito(s0, sd, P, L)
% ItakuraSaito: Itakura-Saito distortion measure
%
% ------
% USAGE: [d, dn] = ItakuraSaito(s0, sd, P, L)
%
%        s0 : original speech signal
%        sd : distorted (processed) speech signal
%        P  : length of the LPC filter [(P-1) is
%             the order of the LPC model]
%        L  : frame size
%        dn : Itakura-Saito measure of each frame
%        d  : average Itakura-Saito measure
%
% ------
% Author:  Dr. Jingdong Chen and Dr. Jacob Benesty
%          Bell Labs, Lucent Technologies
% Copyright 2007(c)
%
if (length(s0) ~= length(sd))
    error('Two signals for comparison have different length.');
end
if (P <= 0)
    error('P must be positive.');
end
if (P >= L)
    error('L must be greater than P.');
end

K = length(s0);

nfrms = fix(K/L);
A0	  = zeros(P,1);
Ad	  = zeros(P,1);

for nf = 1:nfrms
   x0	= zeros(P,1);
   xd	= zeros(P,1);
   R0	= zeros(P,1);
   Rd	= zeros(P,1);
   for i = 1:L
      x0 = [s0((nf-1)*L+i); x0(1:P-1)];
      R0 = R0 + x0(1)*x0;
      %
      xd = [sd((nf-1)*L+i); xd(1:P-1)];
      Rd = Rd + xd(1)*xd;
   end
   R0 = R0/L; Rd = Rd/L; 
   T0 = toeplitz(R0); Td = toeplitz(Rd); 
   T0i = inv(T0); Tdi = inv(Td); 
   uv = zeros(P,1); uv(1) = 1;
   A0 = (T0i*uv)./(uv'*T0i*uv);
   Ad = (Tdi*uv)./(uv'*Tdi*uv);
   
   %Itakura distance
   noma0 = A0'*T0*A0;
   dn(nf) = (Ad'*T0*Ad)/noma0 - 1;
end

dn = dn(:);
d  = mean(dn);
