function diff_noise = gener_diffusenoise(rmic, Lentime, Fs, c, plot_diff_noise)

fs = Fs;
% % Generate sensor signals
P = rmic;
M = length(P(1,:));
% Calculate sensor distances w.r.t. sensor 1
d = zeros(1,M);
for m = 2:M
    d(m) = norm(P(:,m)-P(:,1),2);
end
params.c = c;
params.fs = fs;
params.N_phi = 64;
z = sinf_1D(d,Lentime,params);
diff_noise = z';

if(plot_diff_noise==1)
    % Calculate spatial coherences
    NFFT = 256;     % Number of frequency bins (for analysis)
    w = 2*pi*fs*(0:NFFT/2)/NFFT;
    sc_sim = zeros(M-2,NFFT/2+1);
    sc_theory = zeros(M-2,NFFT/2+1);
    for m = 1:M-1
        [sc,F]=mycohere(z(1,:)',z(m+1,:)',NFFT,fs,hanning(NFFT),0.75*NFFT);
        sc_sim(m,:) = real(sc');
        sc_theory(m,:) = sinc(w*d(m+1)/c/pi);
    end
    
    figure;
    for m = 1:M-1
        color_rgb = rand(3,1);
        plot(F/1000,sc_sim(m,:),'-','Color',color_rgb);hold on;
        plot(F/1000,sc_theory(m,:),'--','Color',color_rgb);hold on;
    end
    xlabel('Frequency [kHz]');
    ylabel('Spatial Coherence');
    title(sprintf('Distance %1.2f m',d(m+1)));
    set(gca,'DataAspectRatio',[1 0.75 1]);
    legend('Simulation','Theory');
    grid on;
end
end