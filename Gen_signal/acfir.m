function [hc, pathnum] = acfir(rs, rmic, lroom, rcwall, npts, highpass)
% ACFIR : Acoustic Channel FIR filter. 
% ------
%    [hc,pathnum] = acfir(rs, rmic, lroom, rcwall, npts, highpass) 
%    calculates the impulse response of an acoustic channel from the 
%    source to the microphone sensor in a rectangluar room by using 
%    the IMAGE method (J.B. Allen & D.A. Berkley, 1979).
% 
% USAGE : [hc,pathnum] = acfir(rs, rmic, lroom, rcwall, npts, highpass)
%	    
%	  rs       : vector radius to the acoustic source in sampling lengths.
%	  rmic     : vector radius to the microphone sensor in sampling lengths.
%	  lroom    : vector of room dimensions in sampling lengths.
%	  rcwall   : matrix of six wall reflection coefficients. (0<=rcwall<=1)
%	  npts     : number of points of the FIR filter to be computed.
%	  highpass : flag of high passing the final impulse response.
%		         = 0, without high-passing;
%		         = 1, with high-passing (default).
%	  hc       : acoustic channel FIR filter impulse response.
%	  pathnum  : Number of paths with the same delay.
%
%	  sampling length = cT, where c is the speed of sound and T is the
%	  		              sampling period.
%
% ------
% Author:  Dr. Jingdong Chen 
% Copyright 2003(c) Bell Labs
%

if (nargin <= 4)
  npts = 1024;
end
if (nargin <= 5)
  highpass = 1;
end

rs = rs(:);
rmic = rmic(:);
lroom = lroom(:);

if (length(rs) ~= 3)
  error('ACFIR: rs demension error.');
end
if (length(rmic) ~= 3)
  error('ACFIR: rmic demension error.');
end
if (length(lroom) ~= 3)
  error('ACFIR: lroom demension error.');
end
if (sum(size(rcwall)==[3 2]) ~= 2)
  error('ACFIR: rcwall demension error.');
end
if (sum(sum(rcwall>=0))~=6 | sum(sum(rcwall<=1))~=6)
  error('ACFIR: wall reflection coefficients have to be in [0,1].');
end
if (npts<=0)
  error('ACFIR: npts<=0 error.');
end

hc = zeros(npts, 1);
pathnum = zeros(npts, 1);
dist = sqrt((rs-rmic).'*(rs-rmic));
if (dist < 0.5)
  hc(1) = 1;
else
  % Find the range of summation.
  nptsx = fix(npts/(2*lroom(1))+1);
  nptsy = fix(npts/(2*lroom(2))+1);
  nptsz = fix(npts/(2*lroom(3))+1);
  for nx = -nptsx:nptsx
    for ly = -nptsy:nptsy
      for mz = -nptsz:nptsz
	imgr = [nx;ly;mz];
	delays = imgdlys(rs, rmic, lroom, imgr);
	index = 0;
	for qx = 0:1
 	  for jy = 0:1
	    for kz = 0:1
	      index = index + 1;
	      intdelay = round(delays(index));
	      if (intdelay <= npts)
	        gid = rcwall(1,1).^abs(nx-qx)*rcwall(1,2).^abs(nx)*...
	              rcwall(2,1).^abs(ly-jy)*rcwall(2,2).^abs(ly)*...
	              rcwall(3,1).^abs(mz-kz)*rcwall(3,2).^abs(mz) ...
		      /(4*pi*delays(index));
	        hc(intdelay) = hc(intdelay) + gid;
		pathnum(intdelay) = pathnum(intdelay) + 1;
	      end
	    end
	  end
	end
      end
    end
  end
end

% Impulse response has been computed. Filter with HPF of 
% 1% of sampling frequency.

% Choice 1: 3rd order Elliptic digital filter
if (highpass ~= 0)
  %[hb, ha] = ellip(4, 0.1, 60, 2/50, 'high');
  %hc = filter(hb, ha, hc);

  % Choice 2: 3rd order IIR LPF from [J.Allen and D.Berkley, 1979].
  %hb = [1.0000; -1.9391; 0.9391];
  %ha = [1.0000; -1.8745; 0.8819];
  %hc = filter(hb, ha, hc);
  
  W = 2.*pi*100.;  T = 1.e-4;
  R1 = exp(-W*T);
  B1 = 2.*R1*cos(W*T); B2 = -R1^2;
  A1 = -(1.+R1); A2 = R1;
  y1 = 0; y2 = 0; y0 = 0;
  for npt=1:npts
    x0 = hc(npt);
    hc(npt)=y0+A1*y1+A2*y2;
    y2 = y1; y1 = y0;
    y0 = B1*y1 + B2*y2 + x0;
  end
end

