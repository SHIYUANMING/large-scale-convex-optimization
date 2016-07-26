clc;clear all;

%cvx_solver sedumi
%cvx_solver sdpt3
%cvx_solver mosek

%cvx_solver scs
cvx_solver_settings('MAX_ITERS', 10^4, 'SCALE', 20);
%cvx_quiet(true)

%%%%%%%%%%%%%
Opt_al=1;
ZF_al=0;
RZF_al=0;
MRT_al=0;

channel_generate=1;
savedata=1;

%%Problem Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=55; K=50; N1=1; Area=10*10^3; N_set=N1*ones(L,1);
P_dB=[-5:10:25];
%P_dB=20;
P=10.^(P_dB./10);
LC=50;              %# loops for channel realizarions

%% Problem Data
params.L=L;   %'L': # RRHs
params.K=K;    %'K': # MUs
params.N_set=N_set;  % set of antennas at all the RRHs
params.delta_set=1*ones(K,1); % set of noise covariance

%% Bi-section algorithm data
params.r_low=0;
%params.r_up=10^3;
params.epsilon=0.01;

%% initialize results
 Maxmin_opt_temp=0; Maxmin_ZFBF_temp=0; Maxmin_MRT_temp=0; Maxmin_RZF_temp=0; 
 
 TransTime_opt_temp=0; SolveTime_opt_temp=0; TotalTime_opt_temp=0;
 TransTime_ZFBF_temp=0; SolveTime_ZFBF_temp=0; TotalTime_ZFBF_temp=0;
 TransTime_MRT_temp=0; SolveTime_MRT_temp=0; TotalTime_MRT_temp=0;
 TransTime_RZF_temp=0; SolveTime_RZF_temp=0; TotalTime_RZF_temp=0;
 
 %% Generate Channel 
 if channel_generate==true
     for ss=1:LC   
         H(:,:,ss)=channel_realization(L, K, N_set, Area);
        % save('H.mat','H');
     end
 end

%load('H.mat');

for loop_channel=1:LC
    
      params.H=H(:,:,loop_channel);

    for loop_power=1:length(P)
        params.P_set=P(loop_power)*ones(L,1);   %set of transmit power constraints at all the RRHs
        
      %% Caculate Upper bound for SINR by Uplink-downlink duality theory
      H_uplink=0;
      for kk=1:K
          H_uplink=H_uplink+params.H(:,kk)*params.H(:,kk)';
      end
      C_sum=log2(det(eye(sum(N_set))+sum(params.P_set)*H_uplink));
      params.r_up=real(2^(C_sum/K))-1;
      
        
    %% Call CVX with Optimal Coordinated Beamforming
    if Opt_al==true
        
        [r, TotalTime, SolveTime]=maxmin_cvx(params);
        
        Maxmin_opt(loop_power)=r;
        TotalTime_opt(loop_power)=TotalTime;
        SolveTime_opt(loop_power)=SolveTime;
        TransTime_opt(loop_power)=TotalTime-SolveTime;
        
    end
    
    %% Call RZF algorithm
    if RZF_al==true
     
         [r, TotalTime, SolveTime]= maxmin_rzf(params);
        
        Maxmin_RZF(loop_power)=r;
        TotalTime_RZF(loop_power)=TotalTime;
        SolveTime_RZF(loop_power)=SolveTime;
        TransTime_RZF(loop_power)=TotalTime-SolveTime;
        
    end
    
    %% Call ZF Beamforming
    if ZF_al==true
        tic;
        [r, TotalTime, SolveTime]= maxmin_zf(params);
        
        Maxmin_ZFBF(loop_power)=r;
        TotalTime_ZFBF(loop_power)=TotalTime;
        SolveTime_ZFBF(loop_power)=SolveTime;
        TransTime_ZFBF(loop_power)=TotalTime-SolveTime;
        
    end
    
    %% Call MRT algorithm
     if MRT_al==true
         tic;
         [r, TotalTime, SolveTime]= maxmin_mrt(params);
         
        Maxmin_MRT(loop_power)=r;
        TotalTime_MRT(loop_power)=TotalTime;
        SolveTime_MRT(loop_power)=SolveTime;
        TransTime_MRT(loop_power)=TotalTime-SolveTime;
        
     end
     
    end
    %% Sumration
    if Opt_al==true  
    Maxmin_opt_temp=Maxmin_opt_temp+Maxmin_opt;
    TransTime_opt_temp=TransTime_opt_temp+TransTime_opt; 
    SolveTime_opt_temp=SolveTime_opt_temp+SolveTime_opt; 
    TotalTime_opt_temp=TotalTime_opt_temp+TotalTime_opt;
    end
    
    if RZF_al==true
    Maxmin_RZF_temp=Maxmin_RZF_temp+Maxmin_RZF;
    TransTime_RZF_temp=TransTime_RZF_temp+TransTime_RZF; 
    SolveTime_RZF_temp=SolveTime_RZF_temp+SolveTime_RZF; 
    TotalTime_RZF_temp=TotalTime_RZF_temp+TotalTime_RZF;
    end
    
    if ZF_al==true
    Maxmin_ZFBF_temp=Maxmin_ZFBF_temp+Maxmin_ZFBF;
    TransTime_ZFBF_temp=TransTime_ZFBF_temp+TransTime_ZFBF; 
    SolveTime_ZFBF_temp=SolveTime_ZFBF_temp+SolveTime_ZFBF; 
    TotalTime_ZFBF_temp=TotalTime_ZFBF_temp+TotalTime_ZFBF;
    end
    
    if MRT_al==true
    Maxmin_MRT_temp=Maxmin_MRT_temp+Maxmin_MRT;
    TransTime_MRT_temp=TransTime_MRT_temp+TransTime_MRT; 
    SolveTime_MRT_temp=SolveTime_MRT_temp+SolveTime_MRT; 
    TotalTime_MRT_temp=TotalTime_MRT_temp+TotalTime_MRT;
    end
end

%% Calculate Average Rate
if Opt_al==true
Maxmin_opt=Maxmin_opt_temp/loop_channel;
TransTime_opt=TransTime_opt_temp/loop_channel;
SolveTime_opt=SolveTime_opt_temp/loop_channel; 
TotalTime_opt=TotalTime_opt_temp/loop_channel;

if savedata==true
save('Maxmin_opt.mat','Maxmin_opt'); save('TransTime_opt.mat','TransTime_opt');
save('SolveTime_opt.mat','SolveTime_opt'); save('TotalTime_opt.mat','TotalTime_opt');
end

end

if RZF_al==true
Maxmin_RZF=Maxmin_RZF_temp/loop_channel;
TransTime_RZF=TransTime_RZF_temp/loop_channel; 
SolveTime_RZF=SolveTime_RZF_temp/loop_channel; 
TotalTime_RZF=TotalTime_RZF_temp/loop_channel; 

if savedata==true
save('Maxmin_RZF.mat','Maxmin_RZF'); save('TransTime_RZF.mat','TransTime_RZF');
save('SolveTime_RZF.mat','SolveTime_RZF'); save('TotalTime_RZF.mat','TotalTime_RZF');
end

end

if ZF_al==true
Maxmin_ZFBF=Maxmin_ZFBF_temp/loop_channel;
TransTime_ZFBF=TransTime_ZFBF_temp/loop_channel; 
SolveTime_ZFBF=SolveTime_ZFBF_temp/loop_channel;
TotalTime_ZFBF=TotalTime_ZFBF_temp/loop_channel;

if savedata==true
save('Maxmin_ZFBF.mat','Maxmin_ZFBF'); save('TransTime_ZFBF.mat','TransTime_ZFBF');
save('SolveTime_ZFBF.mat','SolveTime_ZFBF'); save('TotalTime_ZFBF.mat','TotalTime_ZFBF');
end

end

if MRT_al==true
Maxmin_MRT=Maxmin_MRT_temp/loop_channel;
TransTime_MRT=TransTime_MRT_temp/loop_channel;
SolveTime_MRT=SolveTime_MRT_temp/loop_channel;
TotalTime_MRT=TotalTime_MRT_temp/loop_channel;

if savedata==true
save('Maxmin_MRT.mat','Maxmin_MRT'); save('TransTime_MRT.mat','TransTime_MRT');
save('SolveTime_MRT.mat','SolveTime_MRT'); save('TotalTime_MRT.mat','TotalTime_MRT');
end

end
   
%%
if Opt_al==true
plot(P_dB,log2(1+Maxmin_opt),'b-o','LineWidth',2.0, 'MarkerSize',12); %Optimal beamforming
hold on;
end

if RZF_al==true
plot(P_dB,log2(1+Maxmin_RZF),'b-d','LineWidth',2.0, 'MarkerSize',8); %RZF
hold on;
end

if ZF_al==true
plot(P_dB,log2(1+Maxmin_ZFBF),'g-o','LineWidth',2.0, 'MarkerSize',8); %ZFBF
hold on;
end

if MRT_al==true
plot(P_dB,log2(1+Maxmin_MRT),'y-*','LineWidth',2.0, 'MarkerSize',8); %MRT 
hold on;
end

legend('Optimal Coordinated Beamforming', 'RZF', 'ZFBF', 'MRT');
xlabel('SNR [dB]','fontsize',14,'fontweight','b','fontname','helvetica');
ylabel('Minimum Rate [bps/Hz]','fontsize',14,'fontweight','b','fontname','helvetica');