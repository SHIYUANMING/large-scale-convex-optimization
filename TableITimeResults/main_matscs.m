clc;clear all;
cvx_solver scs

%addpath('/Users/yuanmingshi/Dropbox/Research/Software/scs-master/matlab');
addpath('/home/yuanming/Documents/Software/scs-master/matlab');

%%Problem Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=50; K=50; N1=2; Area=6*10^3; N_set=N1*ones(L,1);
Q=5;
%Q=[0:2:10]';  %QoS in dB
LC=100; %# loops for channel realizarions
prob=zeros(3,1); 
 %"prob(1): stuffing channel realization"
 %"prob(2):stuffing sinr thresholds"
 %"prob(3): stuffing transmit power"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.L=L;   %'L': # RRHs
params.K=K;    %'K': # MUs
params.N_set=N_set;  %set of antennas at all the RRHs
params.delta_set=1*ones(K,1); %set of noise covariance
params.P_set=10^(0)*ones(L,1);   %set of transmit power constraints at all the RRHs


%%%%%%%%%Generate SOCP Structure%%%%%%%%%%%%%
%can be done offline
params.P_set=ones(params.L,1);
params.r_set=ones(params.K,1);  
params.H=ones(sum(params.N_set), params.K)+i*ones(sum(params.N_set), params.K);
[data_socp,cones,dims,data_stuffing] = prob_to_socp(params);

%%%%%%%%%%Parameters for SCS Solver%%%%%%%%%%%%
paramscs.VERBOSE = 0;paramscs.EPS = 1.00e-03; paramscs.ALPHA = 1.8;
paramscs.MAX_ITERS = 10^4; 
paramscs.NORMALIZE = 1; paramscs.SCALE = 20;

%%%%%%initialize results%%%%%%%%%%%%%%%%%%%%%%
SU_counter=0; %count the success of channel realizations
TotalPower_matscs_temp=0;   Time_matscs_temp=0; Timetransfer_matscs_temp=0; Iter_matscs_temp=0;

for ss=1:LC   %generate channel
H(:,:,ss)=channel_realization(L, K, N_set, Area);
end

%load('H.mat');

for lp=1:LC
    
    params.H=H(:,:,lp);
    for lq=1:length(Q)
        params.r_set=10^(Q(lq)/10)*ones(K,1);  %set of SINR thresholds
        
        if lq==1
            prob=[1;1;0];  %stuff for SINR thresholds and channel realizations
        else
            prob=[0;1;0]; %only need to stuff SINR thresholds
        end
%socp data stuffing  
%%
prob=[1;1;1];
tic;
[data_socp, data_stuffing] = matrix_stuffing(params,data_stuffing,prob);
[xx,yy,ss, info] = scs_direct(data_socp,cones,paramscs);
Time_matscs(lq)=toc;

%call SCS solver
%%


if strfind(info.status,'Solved')
    TotalPower_matscs(lq)=norm(xx(dims.n-dims.B+1:dims.n),'fro')^2;
    Timetransfer_matscs(lq)=Time_matscs(lq)-info.solveTime;
    Iter_matscs(lq)=info.iter;
    
else
    TotalPower_matscs(lq)=10^100;
end
    end
 
    if (norm(TotalPower_matscs)<10^100)
    SU_counter=SU_counter+1;
    TotalPower_matscs_temp=TotalPower_matscs_temp+TotalPower_matscs;
    Time_matscs_temp=Time_matscs_temp+Time_matscs;
    Timetransfer_matscs_temp=Timetransfer_matscs_temp+Timetransfer_matscs;
    Iter_matscs_temp=Iter_matscs_temp+Iter_matscs;
    end
end

TotalPower_matscs=TotalPower_matscs_temp/SU_counter
Time_matscs=Time_matscs_temp/SU_counter
Timetransfer_matscs=Timetransfer_matscs_temp/SU_counter
Iter_matscs=Iter_matscs_temp/SU_counter
