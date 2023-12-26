 function gen_impluse_image

%% set up for simulation
% signal parameters
Fs          = 8000*2;   % Sampling frequency in Hz.
impluseLen  = 2048*4;   % length of the impulse responses

% image room for impluse response
c           = 340;                   % Speed of sound in m/s.
sLen        = c/Fs;                  % Sampling length in meter
lroom       = [5; 4; 3];             % Room dimension in meters
rcoeff      = 0.85;                   % wall reflection coefficients

% microphone array positions
M       = 8;
delta   = 0.0425;
rmicy   = 4.0*ones(1,M);
rmicz   = 1.5*ones(1,M);
x_cent  = 3;
x_min   = x_cent - (M-1)*delta/2;
x_max   = x_cent + (M-1)*delta/2;
rmicx   = [x_max : -delta : x_min];
rmic   = [rmicx; rmicy; rmicz]; % microphone location

show_SKAWPE = 1;
if show_SKAWPE
% source location -- 0 degree
rs1 = [5.0, 4.0, 1.5];
% source location -- 180 degree
rs2 = [2.0, 6.0, 1.5];
else
% source location -- 0 degree
rs1 = [5.0, 4.0, 1.5];
% source location -- 180 degree
rs2 = [5.05, 4.0, 1.5];    
end

% show the room layout
figure('position',[0 200 500 300]);
plot(rs1(1),rs1(2),'*b','MarkerSize',10);hold on;
plot(rs2(1),rs2(2),'*r','MarkerSize',10);hold on;
plot(rmic(1,:),rmic(2,:),'ob','MarkerSize',5);hold on;
axis([0 lroom(1) 0 lroom(2)]);box on;
grid on;xlabel('x');ylabel('y');title('setup lavot');


% run imageModel
rcwall = rcoeff * ones(3,2);
if(rcoeff<1e-4)
    highpass = 0; % ideal case without highpass filter
else
    highpass = 1; % reverberation case
end
NoMic = length(rmic(1,:));
H1 = zeros(impluseLen, NoMic);
H2 = zeros(impluseLen, NoMic);
for m = 1 : NoMic
    [hp1, ~] = acfir(rs1/sLen, rmic(:,m)/sLen, lroom/sLen, rcwall, impluseLen, highpass);
    hp1 = hp1/abs(max(hp1));
    H1(:, m) = hp1;
    [hp2, pathcounter] = acfir(rs2/sLen, rmic(:,m)/sLen, lroom/sLen, rcwall, impluseLen, highpass);
    hp2 = hp2/abs(max(hp2));
    H2(:, m) = hp2;
end

save('impulse_response_0_180.mat', 'Fs', 'H1','H2', 'c', 'lroom', 'rcoeff', 'rs1', 'rs2', 'rmic','delta');
% save('impulse_response_0_90.mat', 'Fs', 'H1','H2', 'c', 'lroom', 'rcoeff', 'rs1', 'rs2', 'rmic','delta');
% save('impulse_response_0_135.mat', 'Fs', 'H1','H2', 'c', 'lroom', 'rcoeff', 'rs1', 'rs2', 'rmic','delta');
% save('impulse_response_0_01.mat', 'Fs', 'H1','H2', 'c', 'lroom', 'rcoeff', 'rs1', 'rs2', 'rmic','delta');
  
figure('position',[50, 30, 1140, 550]);
subplot(211);plot(H1);subplot(212);plot(H2);


