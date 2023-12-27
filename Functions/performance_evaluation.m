function PER = performance_evaluation(Para_AWPE, Sig_clean, Sig_observed, Sig_der)

Kl             = Para_AWPE.Kl;
delta_N        = Para_AWPE.delta_N;
alpha_rls      = Para_AWPE.alpha_rls;
Fs             = Para_AWPE.Fs;
cut_kind       = Para_AWPE.cut_kind;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Overall Performance Evaluation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% align amplitude
Sig_clean                 = 0.5*Sig_clean./max(max(Sig_clean));
Sig_observed              = 0.5*Sig_observed./max(max(Sig_observed));
Sig_der                   = 0.5*Sig_der./max(Sig_der);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1) Overall Performance Evaluation

if Para_AWPE.overall_performance
    start_time = 5;
    source          = se_evaluate_score(Sig_clean(start_time*Fs+1 : end),Sig_observed(start_time*Fs+1 : end));
    source_mean     = source(:,1);
    source_median   = source(:,2);
    
    overall_CD_origin       = source_mean(1);
    overall_SRMR_origin     = source_mean(2);
    overall_LLR_origin      = source_mean(3);
    overall_SNR_origin      = source_mean(4);
    overall_PESQ_origin     = source_mean(5);
    
    source2         = se_evaluate_score(Sig_clean(start_time*Fs+1 : end),Sig_der(start_time*Fs+1 : end));
    source_mean2    = source2(:,1);
    source_median2  = source2(:,2);
    
    overall_CD_enhanced     = source_mean2(1);
    overall_SRMR_enhanced   = source_mean2(2);
    overall_LLR_enhanced    = source_mean2(3);
    overall_SNR_enhanced    = source_mean2(4);
    overall_PESQ_enhanced   = source_mean2(5);
    
    PER.overall_CD_origin        = overall_CD_origin;
    PER.overall_LLR_origin       = overall_LLR_origin;
    PER.overall_SNR_origin       = overall_SNR_origin;
    PER.overall_PESQ_origin      = overall_PESQ_origin;
    
    PER.overall_CD_enhanced      = overall_CD_enhanced;
    PER.overall_LLR_enhanced     = overall_LLR_enhanced;
    PER.overall_SNR_enhanced     = overall_SNR_enhanced;
    PER.overall_PESQ_enhanced    = overall_PESQ_enhanced;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % 2) Summary Results %%
    
    Time = datestr(now,31);
    %fprintf('#NR = %d||Kl = %d||delta_N = %d||alpha_la = %.3f||beta_la = %.3f||factor = %.3f#\n',NR,Kl,delta_N,alpha_la,beta_la,robust_factor);
    fprintf('#################################################################\n');
    fprintf('#--------------------- DEREVERBERATION-WPE ---------------------#\n');
    fprintf('#--------------------- Performance Summary ---------------------#\n');
    fprintf('#     Order   = %d    ||  Delay   = %d    || alpha   = %.3f   #\n', Kl, delta_N, alpha_rls);
    fprintf('#---------------The time is: '),fprintf(Time);fprintf('----------------#\n');
    fprintf('#---------------------------------------------------------------#\n');
    fprintf('# Origin CD   = %.3f | Enhanced CD   = %.3f | improved %.3f! #\n', overall_CD_origin, overall_CD_enhanced, overall_CD_origin-overall_CD_enhanced);
    fprintf('# Origin SRMR = %.3f | Enhanced SRMR = %.3f | improved %.3f! #\n', overall_SRMR_origin, overall_SRMR_enhanced, overall_SRMR_enhanced-overall_SRMR_origin);
    fprintf('# Origin LLR  = %.3f | Enhanced LLR  = %.3f | improved %.3f! #\n', overall_LLR_origin, overall_LLR_enhanced, overall_LLR_origin-overall_LLR_enhanced);
    fprintf('# Origin SNR  = %.3f | Enhanced SNR  = %.3f | improved %.3f! #\n', overall_SNR_origin, overall_SNR_enhanced, overall_SNR_enhanced-overall_SNR_origin);
    fprintf('# Origin PESQ = %.3f | Enhanced PESQ = %.3f | improved %.3f! #\n', overall_PESQ_origin, overall_PESQ_enhanced, overall_PESQ_enhanced-overall_PESQ_origin);
    fprintf('#################################################################\n\n');
    
end

if Para_AWPE.segment_performance
    % Doa_change_Info = Fs*Doa_change_Info;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % 3) evaluate performance with small segments
    
    % aviod arrange different speaker's speech in one segment
    Doa_change_Info = [7,17,24,34.5];
    manual_set = 1;
    
    time_vect = cut_segment_eval(cut_kind);
    time_vect = time_vect(1:25,:);
    
    CD_origin       = zeros(length(time_vect(:,1)),1);
    SRMR_origin     = zeros(length(time_vect(:,1)),1);
    LLR_origin      = zeros(length(time_vect(:,1)),1);
    SNR_origin      = zeros(length(time_vect(:,1)),1);
    PESQ_origin     = zeros(length(time_vect(:,1)),1);
    CD_enhanced     = zeros(length(time_vect(:,1)),1);
    SRMR_enhanced   = zeros(length(time_vect(:,1)),1);
    LLR_enhanced    = zeros(length(time_vect(:,1)),1);
    SNR_enhanced    = zeros(length(time_vect(:,1)),1);
    PESQ_enhanced   = zeros(length(time_vect(:,1)),1);
    
    
    for t_index = 1 : length(time_vect(:,1))
        t_start = Fs*time_vect(t_index,1)+1;
        t_end   = Fs*time_vect(t_index,2);
        if t_end > length(Sig_clean)
            t_end   = length(Sig_clean);
        end
        clean_seg                 = Sig_clean(t_start : t_end);
        observed_seg              = Sig_observed(t_start : t_end);
        enhanced_seg              = Sig_der(t_start : t_end);
        
        source          = se_evaluate_score(clean_seg, observed_seg);
        source_mean     = source(:,1);
        source_median   = source(:,2);
        
        CD_origin(t_index)       = source_mean(1);
        SRMR_origin(t_index)     = source_mean(2);
        LLR_origin(t_index)      = source_mean(3);
        SNR_origin(t_index)      = source_mean(4);
        PESQ_origin(t_index)     = source_mean(5);
        
        source2         = se_evaluate_score(clean_seg, enhanced_seg);
        source_mean2    = source2(:,1);
        source_median2  = source2(:,2);
        
        CD_enhanced(t_index)     = source_mean2(1);
        SRMR_enhanced(t_index)   = source_mean2(2);
        LLR_enhanced(t_index)    = source_mean2(3);
        SNR_enhanced(t_index)    = source_mean2(4);
        PESQ_enhanced(t_index)   = source_mean2(5);
    end
    
    PER.time_vect = time_vect;
    PER.CD_origin                   = CD_origin(:);
    PER.LLR_origin                  = LLR_origin(:);
    PER.SNR_origin                  = SNR_origin(:);
    PER.PESQ_origin                 = PESQ_origin(:);
    
    PER.CD_enhanced      = CD_enhanced(:);
    PER.LLR_enhanced     = LLR_enhanced(:);
    PER.SNR_enhanced     = SNR_enhanced(:);
    PER.PESQ_enhanced    = PESQ_enhanced(:);
    
end
