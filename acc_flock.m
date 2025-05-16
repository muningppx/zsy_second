

rate_time=1:0.1:1;
rate_people=1:1:35;
xx_all = cell(length(rate_people),length(rate_time));
%% Set paths

for num=1:length(rate_people)
       o=rate_people(num);
       temp_xx_all_rate=cell(length(rate_time),1);
    for kkkkk=1:length(rate_time)     
        time=rate_time(kkkkk);
        add5 = 'E:/trca/TRCA-SSVEP-master/data/s';
        add6 = num2str(o);
        add7 = '.mat';
        file = strcat(add5, add6, add7);    
        data=load(file);
        
        a=data.a;

         % 假设时间长度�?1�?
   
        fs = 250; % 示例采样频率
   
        eeg = a(:, :, 161:161+time*fs-1, :); 
    
        %% Parameter for analysis (Modify according to your analysis)
        
        % Data length for target identification [s]
        len_gaze_s = 0.5;   
        
        % Visual latency being considered in the analysis [s]
        len_delay_s = 0.13;                  
        
        % The number of sub-bands in filter bank analysis
       
        
        % The number of harmonics in the canonical correlation analysis 
        num_harms = 5;
        
        % 100*(1-alpha_ci): confidence intervals
        alpha_ci = 0.05;                 
        
        %% Fixed parameter (Modify according to the experimental setting)
        
        % Sampling rate [Hz]
        fs = 250;                  
        
        % Duration for gaze shifting [s]
        len_shift_s = 0.5;                  
        
        % List of stimulus frequencies
        list_freqs = [8.6:0.2:15.8 8:0.2:8.4];
                                            
        % The number of stimuli
        num_targs = length(list_freqs);    
        
        % Labels of data
        labels = [1:1:num_targs];         
        
        %% Preparing useful variables (DONT'T need to modify)
        
        % Data length [samples]
        len_gaze_smpl = round(len_gaze_s*fs);           
        
        % Visual latency [samples]
        len_delay_smpl = round(len_delay_s*fs);         
        
        % Selection time [s]
        len_sel_s = len_gaze_s + len_shift_s;
        
        % Confidence interval
        ci = 100*(1-alpha_ci);                  
        is_ensemble=0;
        %% Performing the FBCCA-based SSVEP detection algorithm
        
        fprintf('Results of the FBCCA-based method.\n');
        
        % Preparing data
        


  accs=zeros(6,1);
% Estimate classification performance
for loocv_i = 1:1:6
    a=q{num,kkkkk};
    a=a(loocv_i,:);
    % Training stage 
    traindata = eeg;
    traindata(:, :, :, loocv_i) = [];
    model = train_trca(traindata, fs, a);
    
    % Test stage
    testdata = squeeze(eeg(:, :, :, loocv_i));
    estimated = test_trca(testdata, model, is_ensemble,a);
    estimated_our(loocv_i,:)=estimated;
    % Evaluation 
    is_correct = (estimated==labels);
    accs(loocv_i) = mean(is_correct)*100;
    
    
end % loocv_i
temp_xx_all_rate{kkkkk}=accs;
    end
    xx_all(num,:)=temp_xx_all_rate;
end


function model = train_trca(eeg, fs, x)
target_values_x4 = [2, 3, 4, 5, 6, 7, 8]; % x(3) 要比较的目标�?

% 计算 x(3) 与目标�?�之间的距离
distances_x4 = abs(x(4) - target_values_x4);

% 找到�?近的值的索引
[~, min_index_x4] = min(distances_x4);

% �? x(3) 赋�?�为�?近的�?
x(4) = target_values_x4(min_index_x4);

target_values = [74, 82, 90]; % 要比较的目标�?
% 计算距离
distances = abs(x(5) - target_values);
% 找到�?近的值的索引
[~, min_index] = min(distances);
% 将x(4)赋�?�为�?近的�?
x(5) = target_values(min_index);

[num_targs, num_chans, num_smpls, ~] = size(eeg);
trains = zeros(num_targs, x(4), num_chans, num_smpls);
W = zeros(x(4), num_targs, num_chans);
for targ_i = 1:1:num_targs
    eeg_tmp = squeeze(eeg(targ_i, :, :, :));
    for fb_i = 1:1:x(4)
        eeg_tmp = filterbank(eeg_tmp, fs, fb_i,x(5));
        trains(targ_i,fb_i,:,:) = squeeze(mean(eeg_tmp, 3));
        w_tmp = trca(eeg_tmp);
        W(fb_i, targ_i, :) = w_tmp(:,1);
    end % fb_i
end % targ_i
model = struct('trains', trains, 'W', W,...
     'fs', fs, 'num_targs', num_targs);
end

function W = trca(eeg)

[num_chans, num_smpls, num_trials]  = size(eeg);
S = zeros(num_chans);
for trial_i = 1:1:num_trials-1
    x1 = squeeze(eeg(:,:,trial_i));
    x1 = bsxfun(@minus, x1, mean(x1,2));
    for trial_j = trial_i+1:1:num_trials
        x2 = squeeze(eeg(:,:,trial_j));
        x2 = bsxfun(@minus, x2, mean(x2,2));
        S = S + x1*x2' + x2*x1';
    end % trial_j
end % trial_i
UX = reshape(eeg, num_chans, num_smpls*num_trials);
UX = bsxfun(@minus, UX, mean(UX,2));
Q = UX*UX';
[W,~] = eigs(S, Q);
end

function results = test_trca(eeg, model, is_ensemble,x)


if ~exist('is_ensemble', 'var') || isempty(is_ensemble)
    is_ensemble = 1; end

if ~exist('model', 'var')
    error('Training model based on TRCA is required. See train_trca().'); 
end
target_values_x4 = [2, 3, 4, 5, 6, 7, 8]; % x(3) 要比较的目标�?

% 计算 x(3) 与目标�?�之间的距离
distances_x4 = abs(x(4) - target_values_x4);

% 找到�?近的值的索引
[~, min_index_x4] = min(distances_x4);

% �? x(3) 赋�?�为�?近的�?
x(4) = target_values_x4(min_index_x4);

target_values = [74, 82, 90]; % 要比较的目标�?
% 计算距离
distances = abs(x(5) - target_values);
% 找到�?近的值的索引
[~, min_index] = min(distances);
% 将x(4)赋�?�为�?近的�?
x(5) = target_values(min_index);

fb_coefs = [1:x(4)].^(-x(2))+x(3);

for targ_i = 1:1:model.num_targs
    test_tmp = squeeze(eeg(targ_i, :, :));
    for fb_i = 1:1:x(4)
        testdata = filterbank(test_tmp, model.fs, fb_i,x(5));
        for class_i = 1:1:model.num_targs
            traindata =  squeeze(model.trains(class_i, fb_i, :, :));
            if ~is_ensemble
                w = squeeze(model.W(fb_i, class_i, :));
            else
                w = squeeze(model.W(fb_i, :, :))';
            end
            r_tmp = corrcoef(testdata'*w, traindata'*w);
            r(fb_i,class_i) = r_tmp(1,2);
        end % class_i
    end % fb_i
    rho = fb_coefs*r;
    [~, tau] = max(rho);
    results(targ_i) = tau;
end % targ_i
end

function y = filterbank(eeg, fs, idx_fb,x_5)

if nargin < 2
    error('stats:test_fbcca:LackOfInput', 'Not enough input arguments.'); 
end

if nargin < 3 || isempty(idx_fb)
    warning('stats:filterbank:MissingInput',...
        'Missing filter index. Default value (idx_fb = 1) will be used.'); 
    idx_fb = 1;
elseif idx_fb < 1 || 10 < idx_fb
    error('stats:filterbank:InvalidInput',...
        'The number of sub-bands must be 0 < idx_fb <= 10.'); 
end

[num_chans, ~, num_trials] = size(eeg);
fs=fs/2;

passband = [6, 14, 22, 30, 38, 46, 54, 62, 70, 78];
stopband = [4, 10, 16, 24, 32, 40, 48, 56, 64, 72];
Wp = [passband(idx_fb)/fs, x_5/fs];
Ws = [stopband(idx_fb)/fs, (x_5+10)/fs];
[N, Wn]=cheb1ord(Wp, Ws, 3, 40);
[B, A] = cheby1(N, 0.5, Wn);

y = zeros(size(eeg));
if num_trials == 1
    for ch_i = 1:1:num_chans
        y(ch_i, :) = filtfilt(B, A, eeg(ch_i, :));
    end % ch_i
else
    for trial_i = 1:1:num_trials
        for ch_i = 1:1:num_chans
            y(ch_i, :, trial_i) = filtfilt(B, A, eeg(ch_i, :, trial_i));
        end % trial_i
    end % ch_i
end % if num_trials == 1
end