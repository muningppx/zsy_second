
rate_time=1:0.1:1;
rate_people=1:1:35;
xx_all = cell(length(rate_people),length(rate_time));
parfor num=1:length(rate_people)
       o=rate_people(num);
       temp_xx_all_rate=cell(length(rate_time),1);
    for kkkkk=1:length(rate_time)     
        time=rate_time(kkkkk);
        add5 = 'E:/trca/TRCA-SSVEP-master/data/s';
        add6 = num2str(o);
        add7 = '.mat';
        file = strcat(add5, add6, add7);    
       
        
 
   
        fs = 250; % 绀轰緥閲囨牱棰戠巼
   

    
        %% Parameter for analysis (Modify according to your analysis)
        
        % Data length for target identification [s]
        len_gaze_s = 0.5;   
        
        % Visual latency being considered in the analysis [s]
        len_delay_s = 0.13;                  
        
        % The number of sub-bands in filter bank analysis
        num_fbs = 5;
        
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
        list_freqs = [8:1:15 8.2:1:15.2 8.4:1:15.4 8.6:1:15.6 8.8:1:15.8];
                                            
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
        
        %% Performing the FBCCA-based SSVEP detection algorithm
        
        fprintf('Results of the FBCCA-based method.\n');
        
        % Preparing data
        
        num_blocks=6;
        
        
       
        % Estimate classification performance
        xx=zeros(6,5,5,5);%第三个5_number_第4个5_iter
        x=q{o,kkkkk};
        
% Estimate classification performance
for loocv_i = 1:1:num_blocks
     xx(loocv_i,:,:,:) = ssa( labels,loocv_i,x);
end % loocv_i
    temp_xx_all_rate{kkkkk}=xx;
    disp(kkkkk);
    end
    xx_all(num,:)=temp_xx_all_rate;
end






function acc = test_trca(r,labels,x)
target_values = [2, 3, 4,5,6,7,8]; % 要比较的目标值
% 计算距离
distances = abs(x(3) - target_values);
% 找到最近的值的索引
[~, min_index] = min(distances);

x(3) = target_values(min_index);



target_values = [74, 82, 90]; % 要比较的目标值
% 计算距离
distances = abs(x(4) - target_values);
% 找到最近的值的索引
[~, min_index] = min(distances);
% 将x(4)赋值为最近的值
x(4) = target_values(min_index);

if x(4)==74
    num=1;
elseif x(4)==82
    num=2;
else
    num=3;
end

fb_coefs = [1:x(3)].^(-x(1))+x(2);

for kk =1:5
    for targ_i=1:40
        r_k=squeeze(r(targ_i,1:x(3),:,num,kk));
    
        rho = fb_coefs*r_k;
        [~, tau] = max(rho);
        result(targ_i) = tau;
   

    end
        is_correct = result == labels;
    acc(kk) = mean(is_correct);
end
    acc = mean(acc);
end



function result = ssa(labels,loocv_i,x)
is_ensemble = 0;
if ~exist('is_ensemble', 'var') || isempty(is_ensemble)
    is_ensemble = 1; end


r=squeeze(x(loocv_i,:,:,:,:,:));




for N=10:10:50
    for Max_iter=50:50:250
%N=30;  % 麻雀个数
dim=4;  % 评估函数维度
N_discoverer=0.7*N;  % 发现者个数
N_Followers=0.3*N;   % 追随者个数
SDNumber=0.2*N;   % 警戒者个数
%Max_iter=200;    % 最大迭代次数
ST=0.6;   % 安全阈值

%% 测试函数
f=@(x) 1-test_trca(r,labels,x);

% 边界设置
ub = [1.5, 0.5, 8, 90]; % 上边界
lb = [0, 0, 2, 74];      % 下边界

%% 初始化
x = lb + rand(N,dim) .* (ub - lb);

for i=1:N
    fitness(i) = f(x(i,:));  % 计算麻雀种群的适应度值
end
[A,index]=sort(fitness);
x_best=x(index(1),:);  % 记录所有麻雀走过的位置的最优位置
x_worst=x(index(end),:);  % 记录所有麻雀走过的位置的最差位置
best_fitness=A(1);  % 记录所有麻雀走过的位置的最优值
worst_fitness=A(end);  % 记录所有麻雀走过的位置的最差值
x_best_currently=x(index(1),:);  % 记录当前麻雀种群最优位置
x_worst_currently=x(index(end),:);  % 记录当前麻雀种群最差位置
best_fitness_currently=A(1);   % 记录当前麻雀种群最优值
worst_fitness_currently=A(end);   % 记录当前麻雀种群最差值
x_discoverer=index(1:N_discoverer);   % 发现者位置
x_Followers=index(N_discoverer+1:N);  % 追随者位置
  % 警戒者位置
B=[-1,1];
F=best_fitness;  % 记录每次迭代的麻雀走过的位置的最优值
iter=1;  % 初始化迭代次数

%% 开始迭代更新
while iter<Max_iter
    for i=1:dim
        C(i)=B(round(rand)+1);
    end
    A=C'*inv((C*C'));
    R2=rand;
    % 更新发现者位置
    for i=1:N_discoverer
        for j=1:dim
            if R2<ST
                x(x_discoverer(i),j)=x(x_discoverer(i),j)*exp(-i/(rand*Max_iter));
            else
                 x(x_discoverer(i),j)=x(x_discoverer(i),j)+randn;
            end
        end
        % 边界判断
        ub_flag=x(x_discoverer(i),:)>ub;
        lb_flag=x(x_discoverer(i),:)<lb;
        x(x_discoverer(i),:)=(x(x_discoverer(i),:).*(~(ub_flag+lb_flag)))+ub.*ub_flag+lb.*lb_flag;
    end
     
    % 更新追随者位置
    for i=N_discoverer+1:N
        for j=1:dim
            if i>(N - N_discoverer)/2 + N_discoverer
                x(x_Followers(i-N_discoverer),j)=rand*exp((x_worst_currently(j)-x(x_Followers(i-N_discoverer),j))/i^2);
            else
                x(x_Followers(i-N_discoverer),j)=x_best_currently(j)+abs(x(x_Followers(i-N_discoverer),j)-x_best_currently(j))*A(j);
            end
        end
        % 边界判断
        ub_flag=x(x_Followers(i-N_discoverer),:)>ub;
        lb_flag=x(x_Followers(i-N_discoverer),:)<lb;
        x(x_Followers(i-N_discoverer),:)=(x(x_Followers(i-N_discoverer),:).*(~(ub_flag+lb_flag)))+ub.*ub_flag+lb.*lb_flag;
    end
    % 更新警戒者位置
    Temp = randperm(N);
    SDchooseIndex = Temp(1:SDNumber); 
    for i=1:SDNumber
        for j=1:dim
            if f(x(SDchooseIndex(i),:))~=best_fitness_currently
                x(SDchooseIndex(i),j)=x_best_currently(j)+randn*abs(x(SDchooseIndex(i),j)-x_best_currently(j));
            else
                x(SDchooseIndex(i),j)=x(SDchooseIndex(i),j)+B(round(rand)+1)*(abs(x(SDchooseIndex(i),j)-x_worst_currently(j)))/abs(f(x(SDchooseIndex(i),:))-worst_fitness_currently)+1;
            end
        end
        % 边界判断
        ub_flag=x(SDchooseIndex(i),:)>ub;
        lb_flag=x(SDchooseIndex(i),:)<lb;
        x(SDchooseIndex(i),:)=(x(SDchooseIndex(i),:).*(~(ub_flag+lb_flag)))+ub.*ub_flag+lb.*lb_flag;
    end
    
    for i=1:N
        fitness(i)=f(x(i,:));   % 计算适应度
    end
    [E,index]=sort(fitness);
    if f(x(index(1),:))<best_fitness    % 更新所有麻雀走过的位置的最优位置和最优值
        best_fitness=f(x(index(1),:));
        x_best=x(index(1),:);
    end
    if f(x(index(end),:))>worst_fitness   % 更新所有麻雀走过的位置的最差位置和最差值
        worst_fitness= f(x(index(end),:));
        x_worst=x(index(end),:);
    end
    x_best_currently=x(index(1),:);   % 更新当前麻雀种群的最优位置
    x_worst_currently=x(index(end),:);  % 更新当前麻雀种群的最差位置
    best_fitness_currently=E(1);   % 更新当前麻雀的种群的最优值
    worst_fitness_currently=E(end);  % 更新当前麻雀的种群的最差值
    x_discoverer=index(1:N_discoverer);  % 重新选择种群中的发现者
    x_Followers=index(N_discoverer+1:N_discoverer+N_Followers);  % 重新选择种群中的追随者
   % 重新选择种群中的警戒者
    F(iter)=best_fitness_currently;  % 记录每次迭代的最优值
    iter=iter+1;   % 迭代次数加1
end
result(:,N/10,Max_iter/50)=[F(end) x_best];
    end
end
% display(['最优值是:',num2str(F(end)),'最优麻雀位置:',num2str(x_best)]);
% figure(1);
% plot(F);
% xlabel('迭代次数'),ylabel('适应度值');
end
