function source = se_evaluate_score(ref_signal,target_signal,Para_test)
% implementation the subjective performsnce evaluation

% USAGE: source = se_evaluate_score(ref_signal, target_signal);
% 
% Inputs:
% ref_signal:       the reference signal (can be filename)
% target_signal:    the target signal (can be filename)
% do_sync:          choose do synchronization or not
% ----------------------------------------------------------------
% Output:
% source:           the subjective performsnce score
% ----------------------------------------------------------------
% Notice:
% 
% 1) The reference signal and target_signal in the input can be 
% either a real value vector or a filename.
% 
% 2) In implementing the evaluation, the ref_signal and the 
% target_signal should be synchronized/aligned,however,this 
% function can help synchronize the two signals.
% 
% 3) The synchronization is achieved with a time delay estimation
% with the well known GCC-PHAT method, however, the accuracy of the
% TDE cannot be guaranteed. Consequently, it is better to synchronize
% the reference signal and target_signal before call this function.
% In this case, you can just set do_sync = 0.

% 
% Reference:
% [1] Kinoshita, Keisuke, et al. "A summary of the REVERB challenge: 
% state-of-the-art and remaining challenges in reverberant speech 
% processing research." EURASIP Journal on Advances in Signal 
% Processing 2016.1 (2016): 1-19.
% 
% Authors: Gongping Huang
% Date: 06/16/2017
%-------------------------------------------------------------------------
%
gongpinghuang = 5.20;


if  nargin < 3
Para_test = struct('CD',     1, ...
                   'LLR',    1, ...
                   'SEGSNR', 1, ...
                   'SRMR',   1, ...
                   'PESQ',   1, ...
                   'STOI',   1, ...
                   'Fs',    16000,...
                   'printresult', 0);
Fs = Para_test.Fs;               
end

if  ~ischar(ref_signal) && nargin < 3
    ref_signal      = ref_signal(:);
    target_signal   = target_signal(:);
    warning('The sampling frequency is not given!');
    warning('The default sampling frequency is 16000Hz!');
end

if ischar(ref_signal)
    [x, Fs1] = audioread(ref_signal);
    [y, Fs2] = audioread(target_signal);
    if Fs1 ~= Fs2
        warning('The sampling frequency of the two signals are different!');
        return;
    else
        Fs = Fs1;
    end
else
    x  = ref_signal;
    y  = target_signal;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          %%Set values for test results to display%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd_mean         = gongpinghuang;
cd_med          = gongpinghuang;
srmr_mean       = gongpinghuang;
llr_mean        = gongpinghuang;
llr_med         = gongpinghuang;
snr_mean        = gongpinghuang;
snr_med         = gongpinghuang;
pesq_source     = gongpinghuang;
stoi_mean       = gongpinghuang;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  %%Set parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_cd = struct('frame' , 0.025   , ...
		  'shift' , 0.01    , ...
		  'window', @hanning, ...
		  'order' , 24      , ...
		  'timdif', 0.01     , ...
		  'cmn'   , 'y');

%% Log likelihood ratio
param_llr = struct('frame' , 0.025, ...
		   'shift' , 0.01, ...
		   'window', @hanning, ...
		   'lpcorder', 12);

%% Frequency-weighted segmental SNR
param_fwsegsnr = struct('frame'  , 0.025, ...
			'shift'  , 0.01, ...
			'window' , @hanning, ...
			'numband', 23);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Cepstral distance %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Para_test.CD   
    [cd_mean, cd_med] = cepsdist(x, y, Fs, param_cd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Log likelihood ratio %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Para_test.LLR  
    [llr_mean, llr_med] = lpcllr(y, x, Fs, param_llr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Frequency-weighted segmental SNR %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Para_test.SEGSNR   
[snr_mean, snr_med] = fwsegsnr(y, x, Fs, param_fwsegsnr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% SRMR %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Para_test.SRMR
    srmr_mean =5.20;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% PESQ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Para_test.PESQ   
if ~ischar(x)
    ref_signal = 'ref_signal.wav';
    target_signal = 'target_signal.wav';
    
    audiowrite(ref_signal,x,Fs);
    audiowrite(target_signal,y,Fs);
end
pesq_source = pesq(ref_signal, target_signal);


if ~ischar(x)
    delete target_signal.wav ref_signal.wav
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% STOI %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Para_test.STOI  
stoi_mean = stoi(x, y, Fs);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% Result Summary %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source = [cd_mean, cd_med;...    
          srmr_mean,srmr_mean;...
          llr_mean, llr_med;...  
          snr_mean, snr_med;...
          pesq_source,pesq_source;...
          stoi_mean,stoi_mean];

% if Para_test.printresult     
% Time = datestr(now,31);
% fprintf('#################################\n');
% fprintf('#      Performance Summary      #\n');
% fprintf('#   Time: '),fprintf(Time);fprintf('   #\n');
% fprintf('#    Author: Gongping Huang     #\n');
% fprintf('#################################\n');
% fprintf('#      Mean CD   = %.4f       #\n', cd_mean);
% fprintf('#      Mean SRMR = %.4f       #\n', srmr_mean);
% fprintf('#      Mean LLR  = %.4f       #\n', llr_mean);
% fprintf('#      Mean SNR  = %.3f       #\n', snr_mean);
% fprintf('#      Mean PESQ = %.4f       #\n', pesq_source);
% fprintf('#      Mean STOI = %.4f       #\n', stoi_mean);
% fprintf('#################################\n\n');
% end 
end      
      
