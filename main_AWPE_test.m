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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1) load the noisy signal   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Parameters for dereverberation (DR)
M                       = 8;
L_vect = 4;
alpha_rls = 0.995;
%lamuda_min_vect=5;
%perfermance = zeros(length(L_vect),4);
perfermance = zeros(length(L_vect),4);

for L_index=1:length(L_vect)
    Kl                      = L_vect(L_index);
    %Kl                      = 8;
    
    %alpha_rls               = alpha_rls_vect(alpha_index); % forgetting factor for RLS adaptive
    Para_AWPE.Kl            = Kl;
    Para_AWPE.alpha_rls     = alpha_rls;
    % minimum value of the variance on the denominator
    Para_AWPE.lamuda_min    = 10^(-5); 
    
    load('cleansignal_16k.mat');
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
        delta_N           = 20;        % the delay in frame
        PM_data.x_max     = 40;
        PM_data.CD_min    = 0; PM_data.CD_max = 5;
        PM_data.LLR_min   = 0; PM_data.LLR_max = 0.8;
        PM_data.SNR_min   = 0; PM_data.SNR_max = 25;
        PM_data.PESQ_min  = 0; PM_data.PESQ_max = 4.5;
    end
    Para_AWPE.delta_N        = delta_N;

    Para_AWPE.ct = 0.5;

    % % call the  AWPE dereverbeation filter

    PM_data.CD_matrix   = [];
    PM_data.LLR_matrix  = [];
    PM_data.SNR_matrix  = [];
    PM_data.PESQ_matrix = [];

    DE_method  = 'AWPE';
    [Sig_der, delay] = derev_awpe_RLS_DR(y,Para_AWPE);
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
    
    perfermance(L_index,1) = PER.overall_CD_origin-PER.overall_CD_enhanced;
    perfermance(L_index,2) = PER.overall_LLR_origin-PER.overall_LLR_enhanced;
    perfermance(L_index,3) = PER.overall_SNR_enhanced-PER.overall_SNR_origin;
    perfermance(L_index,4) = PER.overall_PESQ_enhanced-PER.overall_PESQ_origin;

end

fig=figure('position',[0 0 1000 600]);
linewd          = 1.5;
FontSize        = 12;
subplot(221);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(L_vect, perfermance(:,1),'-sb','LineWidth',linewd);grid on;
xlabel('L','fontsize',FontSize,'interpreter', 'latex')
ylabel('$\Delta$CD','fontsize',FontSize,'interpreter', 'latex')
%text(0.5*PM_data.x_max,PM_data.CD_max+0.5,'(a)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
%axis([0 length(time_vect) PM_data.CD_min PM_data.CD_max]);box on;grid on;

subplot(222);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(L_vect, perfermance(:,2),'-sb','LineWidth',linewd);grid on;
xlabel('L','fontsize',FontSize,'interpreter', 'latex')
ylabel('$\Delta$LLR','fontsize',FontSize,'interpreter', 'latex')
%text(0.5*PM_data.x_max,PM_data.LLR_max+0.1,'(b)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
%axis([0 length(time_vect) PM_data.LLR_min PM_data.LLR_max]);box on;grid on;

subplot(223);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(L_vect, perfermance(:,3),'-sb','LineWidth',linewd);grid on;
xlabel('L','fontsize',FontSize,'interpreter', 'latex')
ylabel('$\Delta$FWSNR(dB)','fontsize',FontSize,'interpreter', 'latex')
%text(0.5*PM_data.x_max,PM_data.SNR_max+2,'(c)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
%axis([0 length(time_vect) PM_data.SNR_min PM_data.SNR_max]);box on;grid on;

subplot(224);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(L_vect, perfermance(:,4),'-sb','LineWidth',linewd);grid on;
xlabel('L','fontsize',FontSize,'interpreter', 'latex')
ylabel('$\Delta$PESQ','fontsize',FontSize,'interpreter', 'latex')
%text(0.5*PM_data.x_max,PM_data.PESQ_max+0.5,'(d)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
%axis([0 length(time_vect) PM_data.PESQ_min PM_data.PESQ_max]);box on;grid on;


%{
P = length(PM_data.CD_matrix(1,:));

CD_origin   = repmat(PM_data.CD_origin, 1, P);
LLR_origin  = repmat(PM_data.LLR_origin, 1, P);
SNR_origin  = repmat(PM_data.SNR_origin, 1, P);
PESQ_origin  = repmat(PM_data.PESQ_origin, 1, P);

CD_matrix   = CD_origin-PM_data.CD_matrix;
LLR_matrix  = LLR_origin-PM_data.LLR_matrix;
SNR_matrix  = PM_data.SNR_matrix-SNR_origin;
PESQ_matrix = PM_data.PESQ_matrix-PESQ_origin;

time_vect = time_vect(:,2);

fig=figure('position',[0 0 1000 600]);
linewd          = 1.5;
FontSize        = 12;
subplot(221);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(time_vect, CD_matrix,'-ob','LineWidth',linewd);grid on;
xlabel('Time (s)','fontsize',FontSize,'interpreter', 'latex')
ylabel('$\Delta$CD','fontsize',FontSize,'interpreter', 'latex')
%text(0.5*PM_data.x_max,PM_data.CD_max+0.5,'(a)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 length(time_vect) PM_data.CD_min PM_data.CD_max]);box on;grid on;

subplot(222);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(time_vect, LLR_matrix,'-ob','LineWidth',linewd);grid on;
xlabel('Time (s)','fontsize',FontSize,'interpreter', 'latex')
ylabel('$\Delta$LLR','fontsize',FontSize,'interpreter', 'latex')
%text(0.5*PM_data.x_max,PM_data.LLR_max+0.1,'(b)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 length(time_vect) PM_data.LLR_min PM_data.LLR_max]);box on;grid on;

subplot(223);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(time_vect, SNR_matrix,'-ob','LineWidth',linewd);grid on;
xlabel('Time (s)','fontsize',FontSize,'interpreter', 'latex')
ylabel('$\Delta$FWSNR(dB)','fontsize',FontSize,'interpreter', 'latex')
%text(0.5*PM_data.x_max,PM_data.SNR_max+2,'(c)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 length(time_vect) PM_data.SNR_min PM_data.SNR_max]);box on;grid on;

subplot(224);
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(time_vect, PESQ_matrix,'-ob','LineWidth',linewd);grid on;
xlabel('Time (s)','fontsize',FontSize,'interpreter', 'latex')
ylabel('$\Delta$PESQ','fontsize',FontSize,'interpreter', 'latex')
%text(0.5*PM_data.x_max,PM_data.PESQ_max+0.5,'(d)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 length(time_vect) PM_data.PESQ_min PM_data.PESQ_max]);box on;grid on;
%}