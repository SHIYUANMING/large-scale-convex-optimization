%This Matlab script can be used to generate the results by "Matrix Stuffing+SCS" in Table I in the
%paper:

%Yuanming Shi, Jun Zhang, Brendan O'Donoghue, and Khaled B. Letaief, "Large-scale convex 
%optimization for dense wireless cooperative networks," IEEE Trans. Signal Process., to appear, 2015.

%The implementation of this algorithm needs to download the solver "scs" 
% https://github.com/cvxgrp/scs 
%based on the paper

%Brendan O'Donoghue, Eric Chu, Neal Parikh, and Stephen Boyd, "Conic optimization 
%via operator splitting and  homogeneous self-dual embedding," arXiv preprint arXiv:1312.3039, 2013.

%This is version 1.0 (Last edited: 2015-06-01)

%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.

%Note: the channels are generated randomly, thus the results
%will not be exactly the same as in the paper.


clc;clear all;

%add the path for the solver "scs"
addpath('/scs-master/matlab');
addpath(genpath(pwd));

%%Problem Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=50; K=50; N1=2; Area=6*10^3; N_set=N1*ones(L,1); %'L': # RAU; 'K': # MUs; 'N1: # antennas in each RAU'
Q=5; %QoS in dB
LC=1; %# loops for channel realizarions
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
[data_socp,cones,dims,data_stuffing] = prob_to_socp(params);  % generate and store standard form

%%%%%%%%%%Parameters for SCS Solver%%%%%%%%%%%%
paramscs.VERBOSE = 0;paramscs.EPS = 1.00e-03; paramscs.ALPHA = 1.8;
paramscs.MAX_ITERS = 10^4; 
paramscs.NORMALIZE = 1; paramscs.SCALE = 20;

%%%%%%initialize results%%%%%%%%%%%%%%%%%%%%%%
SU_counter=0; %count the success of channel realizations
TotalPower_matscs_temp=0;  %objective 
Time_matscs_temp=0; % total time
Timetransfer_matscs_temp=0; %modeling time

for ss=1:LC   %generate channel
H(:,:,ss)=channel_realization(L, K, N_set, Area);
end

for lp=1:LC
    
    params.H=H(:,:,lp);
    for lq=1:length(Q)
        params.r_set=10^(Q(lq)/10)*ones(K,1);  %set of SINR thresholds
        
        if lq==1
            prob=[1;1;0];  %stuff for SINR thresholds and channel realizations
        else
            prob=[0;1;0]; %only need to stuff SINR thresholds
        end
        
  
%%solve problem instances 
prob=[1;1;1];
tic;
[data_socp, data_stuffing] = matrix_stuffing(params,data_stuffing,prob); %standard form data stuffing
[xx,yy,ss, info] = scs_direct(data_socp,cones,paramscs); %call SCS solver
Time_matscs(lq)=toc;


%%Prepare outputs

if strfind(info.status,'Solved')
    TotalPower_matscs(lq)=norm(xx(dims.n-dims.B+1:dims.n),'fro')^2;
    Timetransfer_matscs(lq)=Time_matscs(lq)-info.solveTime;
else
    TotalPower_matscs(lq)=10^100;
end
    end
 
    if (norm(TotalPower_matscs)<10^100)
    SU_counter=SU_counter+1;
    TotalPower_matscs_temp=TotalPower_matscs_temp+TotalPower_matscs;
    Time_matscs_temp=Time_matscs_temp+Time_matscs;
    Timetransfer_matscs_temp=Timetransfer_matscs_temp+Timetransfer_matscs;
    end
end

TotalPower_matscs=TotalPower_matscs_temp/SU_counter % objective
Time_matscs=Time_matscs_temp/SU_counter % total time
Timetransfer_matscs=Timetransfer_matscs_temp/SU_counter % modeling time
