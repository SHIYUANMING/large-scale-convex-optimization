clc;clear all;
cvx_solver sedumi
%cvx_solver sdpt3
%cvx_solver mosek

 %cvx_solver scs
 %cvx_solver_settings('NORMALIZE', 1, 'SCALE',20, 'MAX_ITERS', 10^4)
 %cvx_quiet(true)

%addpath('/Users/yuanmingshi/scs/');

%%Problem Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=50; K=50; N1=2; Area=6*10^3; N_set=N1*ones(L,1);
Q=5;
%Q=[0:2:10]';  %QoS in dB
LC=100;  % # loops for channel realizarions
prob=zeros(3,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.L=L;   %'L': # RRHs
params.K=K;    %'K': # MUs
params.N_set=N_set;  %set of antennas at all the RRHs
params.delta_set=1*ones(K,1); %set of noise covariance
params.P_set=10^(0)*ones(L,1);   %set of transmit power constraints at all the RRHs

%%%%%%initialize results%%%%%%%%%%%%%%%%%%%%%%
SU_counter=0; %count the success of channel realizations
TotalPower_matscs_temp=0;   TotalPower_cvx_temp=0;
Time_matscs_temp=0; Time_cvx_temp=0;
Iter_cvx_temp=0; 
TotalTime_cvx_temp=0;
SolveTime_cvx_temp=0;
TransTime_cvx_temp=0;

for ss=1:LC   %generate channel
H(:,:,ss)=channel_realization(L, K, N_set, Area);
end
%   load('H.mat');

for lp=1:LC
    params.H=H(:,:,lp);
    
    for lq=1:length(Q)
        params.r_set=10^(Q(lq)/10)*ones(K,1);  %set of SINR thresholds   
%%%%%%%%%%CVX%%%%%%%%%%%%%%%%%%%%%%%
%%Use Software CVX
tic;
%[Wsolution_cvx, feasible_cvx, solving_time, cvx_slvitr] = powermin_cvx(params);
[Wsolution_cvx, feasible_cvx] = powermin_cvx(params);
TotalTime_cvx(lq)=toc;
%SolveTime_cvx(lq)=solving_time;
%TransTime_cvx(lq)=TotalTime_cvx(lq)-SolveTime_cvx(lq);
%Iter_cvx(lq)=cvx_slvitr;

if feasible_cvx==1
    TotalPower_cvx(lq)=norm(Wsolution_cvx,'fro')^2;
else
    TotalPower_cvx(lq)=10^100;
end
    end
 
    if (norm(TotalPower_cvx)<10^100)
    SU_counter=SU_counter+1;
    TotalPower_cvx_temp=TotalPower_cvx_temp+TotalPower_cvx;
    TotalTime_cvx_temp=TotalTime_cvx_temp+TotalTime_cvx;
    %SolveTime_cvx_temp=SolveTime_cvx_temp+SolveTime_cvx;
    %TransTime_cvx_temp=TransTime_cvx_temp+TransTime_cvx;
    %Iter_cvx_temp=Iter_cvx_temp+Iter_cvx;
    end
end

TotalPower_cvx=TotalPower_cvx_temp/SU_counter
TotalTime_cvx=TotalTime_cvx_temp/SU_counter
%SolveTime_cvx=SolveTime_cvx_temp/SU_counter
%TransTime_cvx=TransTime_cvx_temp/SU_counter
%Iter_cvx=Iter_cvx_temp/SU_counter



