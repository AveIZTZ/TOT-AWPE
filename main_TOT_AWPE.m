% This is to run the recursive (forgetting factor) based adaptive
% weighted prediction error (AWPE) algorithm for dereverberation
% (speech enhancement).
%
%
% Authors: Yujie Zhu
% Date: 08/08/2023
%-------------------------------------------------------------------------

clear;clc;
close all;

% all path in the current folder are included.
warning('OFF')
addpath(genpath(pwd))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1) often change parameters   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_16k_signal          = 1;
overall_performance     = 1;
segment_performance     = 1;

% % choose to use direct path signal or early reflections as reference signal
ref_as_early_refelection= 1; % 1: early reflection; 0) direct path
% length of the start seconds for initalization
% Para_AWPE.init_time     = 5.0; % with italization
Para_AWPE.init_time   = 0.001; % without italization

% Parameters for dereverberation (DR)
M                       = 8;
Kl                      = 16;
alpha_rls               = 0.995; % forgetting factor for RLS adaptive
Para_AWPE.Kl            = Kl;
Para_AWPE.alpha_rls     = alpha_rls;
% minimum value of the variance on the denominator
Para_AWPE.lamuda_min    = 1e-3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1) load the noisy signal   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load('cleansignal_0_01_16k.mat');
load('cleansignal_16k.mat');
%load('cleansignal_0_135_16k_change.mat');
%load('Measured_Signal_t610_snr20_0_30_16k_change.mat');
%load('Measured_Signal_t610_0_30_16k_change.mat');
%time_range = [0*Fs+1:40.1*Fs];
time_range = [0*Fs+1:25*Fs];
Para_AWPE.cut_kind = 2;

if ref_as_early_refelection
    x      = Xearly(time_range,1); % only consider the first channel
else
    x      = Sclean(time_range,1);
end

y          = Xclean(time_range,1:M);
y          = 0.5*y./max(max(y));
M          = length(y(1,:));

% clear some data and parameters
clear Xclean Sclean Xearly time_range

Para_AWPE.Fs          = Fs;
% % choose to use direct path signal or early reflections as reference signal
if ref_as_early_refelection
    delta_N           = 10;        % the delay in frame
    PM_data.x_max     = 40;
    PM_data.CD_min    = 0; PM_data.CD_max = 5;
    PM_data.LLR_min   = 0; PM_data.LLR_max = 0.8;
    PM_data.SNR_min   = 0; PM_data.SNR_max = 25;
    PM_data.PESQ_min  = 0; PM_data.PESQ_max = 4.5;
else
    delta_N         = 2;        % the delay in frame
    PM_data.x_max     = 40;
    PM_data.CD_min    = 0; PM_data.CD_max = 5;
    PM_data.LLR_min   = 0; PM_data.LLR_max = 0.8;
    PM_data.SNR_min   = 0; PM_data.SNR_max = 25;
    PM_data.PESQ_min  = 0; PM_data.PESQ_max = 4.5;
end
Para_AWPE.delta_N        = delta_N;

P_vect = [1,2,3,4];

Para_AWPE.ct = 0.5;

% % call the  AWPE dereverbeation filter

PM_data.CD_matrix   = [];
PM_data.LLR_matrix  = [];
PM_data.SNR_matrix  = [];
PM_data.PESQ_matrix = [];

DE_method  = 'AWPE';
[Sig_der, delay] = derev_awpe_RLS(y,Para_AWPE);
delay            = delay+1;
Sig_der          = Sig_der(delay : end-delay);
% length
Sig_clean        = x(1:length(Sig_der),1);
Sig_observed     = y(1:length(Sig_der),1);
    
%clear x y DE_method  Doa_Info 
% % Implement Performance evaluation
Para_AWPE.overall_performance = overall_performance;
Para_AWPE.segment_performance = segment_performance;
PER = performance_evaluation(Para_AWPE, Sig_clean, Sig_observed, Sig_der);
    
time_vect = PER.time_vect;                                              
PM_data.CD_origin                   = PER.CD_origin(:);
PM_data.LLR_origin                  = PER.LLR_origin(:);
PM_data.SNR_origin                  = PER.SNR_origin(:);
PM_data.PESQ_origin                 = PER.PESQ_origin(:);
    
PM_data.CD_matrix(:,1)      = PER.CD_enhanced(:);
PM_data.LLR_matrix(:,1)     = PER.LLR_enhanced(:);
PM_data.SNR_matrix(:,1)     = PER.SNR_enhanced(:);
PM_data.PESQ_matrix(:,1)    = PER.PESQ_enhanced(:);

for p_index = 1:length(P_vect)
    Para_AWPE.P       = P_vect(p_index);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % 2) Implement AWPE Dereverberation %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DE_method  = 'AWPE-TOT';
    Para_AWPE.L11      = 8;
    Para_AWPE.L12      = 8;
    Para_AWPE.L2      = 2;
    [Sig_der, delay] = derev_awpe_RLS_TOT(y, Para_AWPE);   
    
    delay            = delay+1;
    Sig_der          = Sig_der(delay : end-delay);
    % length
    Sig_clean        = x(1:length(Sig_der),1);
    Sig_observed     = y(1:length(Sig_der),1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % 2)      Dereverberation Performance Evaluation            %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % % Implement Performance evaluation
    Para_AWPE.overall_performance = overall_performance;
    Para_AWPE.segment_performance = segment_performance;
    PER = performance_evaluation(Para_AWPE, Sig_clean, Sig_observed, Sig_der);
    
    time_vect = PER.time_vect;
    %PM_data.CD_origin                   = PER.CD_origin(:);
    %PM_data.LLR_origin                  = PER.LLR_origin(:);
    %PM_data.SNR_origin                  = PER.SNR_origin(:);
    %PM_data.PESQ_origin                 = PER.PESQ_origin(:);
    
    PM_data.CD_matrix(:,p_index+1)      = PER.CD_enhanced(:);
    PM_data.LLR_matrix(:,p_index+1)     = PER.LLR_enhanced(:);
    PM_data.SNR_matrix(:,p_index+1)     = PER.SNR_enhanced(:);
    PM_data.PESQ_matrix(:,p_index+1)    = PER.PESQ_enhanced(:);

end

% save data_KAWPE_P_2Mic PM_data time_vect

% save data_AWPE_switch_simu PM_data time_vect Sig_clean Sig_observed Sig_der
Fun_plot_figure_method_TOT(time_vect(:,2),PM_data)
save('PM_data.mat', 'PM_data','time_vect');

Fun_plot_spectrum_AWPE(Sig_clean, Sig_observed, Sig_der, Fs);
