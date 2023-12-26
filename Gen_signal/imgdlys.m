function delays=imgdlys(rs, rmic, lroom, imgr)
% IMGDLYS : Delays of eight images of a point in a rectangle room.
% ------
%     delays=imgdlys(rs, rmic, lroom, imgr) computes delays of 
%     eight images of a point in a rectangle room.
%     
% USAGE : delays=imgdlys(rs, rmic, lroom, imgr)
%
%         rs     : vector radius to the acoustic source.
%         rmic   : vector radius to the microphone sensor.
%         lroom  : vector of room dimensions.
%	      imgr   : vector of mean image number [nx;ly;mz].
%	      delays : vector of delays in sampling periods (cT) of
%		           eight images of the given point.
%
% ------
% Author:  Dr. Jingdong Chen
% Copyright 2003(c) Bell Labs
%
index = 0;
for px = -1:2:1
  for py = -1:2:1
    for pz = -1:2:1
      index = index + 1;
      Rp(:,index) = rs(:) + [px;py;pz].*rmic(:);
    end
  end
end
Rr = 2*imgr(:).*lroom(:);
delays = sqrt(sum((Rp+Rr*ones(1,8)).^2));
delays = delays(:);

