%This Matlab script can be used to generate the results by "Matrix
%Stuffing+SCS" in Fig.5 for max-min rate optimization in the paper:

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

%%add the path for the solver "scs"
addpath('/scs-master/matlab');
addpath(genpath(pwd));

%%Problem Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=55; K=50; N1=1; Area=10*10^3; N_set=N1*ones(L,1);  %'L': # RAU; 'K': # MUs; 'N1: # antennas in each RAU'
P_dB=[-5:10:25];  %QoS in dB
P=10.^(P_dB./10);
LC=1;              %# loops for channel realizarions

prob=zeros(3,1);  %"prob(1): stuffing channel realization"  %"prob(2):stuffing sinr thresholds" %"prob(3): stuffing transmit power"

%% Problem Data
params.L=L;   %'L': # RRHs
params.K=K;    %'K': # MUs
params.N_set=N_set;  %set of antennas at all the RRHs
params.delta_set=1*ones(K,1); %set of noise covariance

%% Generate SOCP Structure (Can be done offline)
params.P_set=ones(params.L,1);
params.r_set=ones(params.K,1);  
params.H=ones(sum(params.N_set), params.K)+i*ones(sum(params.N_set), params.K);
[data_socp,cones,dims,data_stuffing] = prob_to_socp(params); % generate and store standard form

%% Bi-section algorithm data
params.r_low=0; % rate lower bound
params.epsilon=0.01; 

%% initialize results
 Maxmin_matscs_temp=0;  
 TransTime_matscs_temp=0; SolveTime_matscs_temp=0; TotalTime_matscs_temp=0;
 
 %% Generate Channel 
for ss=1:LC   
H(:,:,ss)=channel_realization(L, K, N_set, Area);
end

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
      params.r_up=real(2^(C_sum/K))-1;  % rate upper bound
      
    %% Matrix Stuffing with SCS for Optimal Coordinated Beamforming
    
     [r, TransTime, SolveTime]=maxmin_matscs(params, data_stuffing, cones); 
     
     Maxmin_matscs(loop_power)=r; % rate
     TransTime_matscs(loop_power)=TransTime; % modeling time
     SolveTime_matscs(loop_power)=SolveTime; % solving time
     TotalTime_matscs(loop_power)=TransTime+SolveTime; % total time
     
    end
    %%
    Maxmin_matscs_temp=Maxmin_matscs_temp+Maxmin_matscs;
    TransTime_matscs_temp=TransTime_matscs_temp+TransTime_matscs; 
    SolveTime_matscs_temp=SolveTime_matscs_temp+SolveTime_matscs; 
    TotalTime_matscs_temp=TotalTime_matscs_temp+TotalTime_matscs;
end
   
%%Prepare outputs 
   Maxmin_matscs=Maxmin_matscs_temp/loop_channel; % rate
   TransTime_matscs=TransTime_matscs_temp/loop_channel; % modeling time
   SolveTime_matscs=SolveTime_matscs_temp/loop_channel; % solving time
   TotalTime_matscs=TotalTime_matscs_temp/loop_channel; % total time

%%
plot(P_dB,log2(1+Maxmin_matscs),'g-d','LineWidth',2.0, 'MarkerSize',15); %Operator splitting method
hold on;

legend('Optimal Coordinated Beamforming');
xlabel('SNR [dB]','fontsize',14,'fontweight','b','fontname','helvetica');
ylabel('Minimum Rate [bps/Hz]','fontsize',14,'fontweight','b','fontname','helvetica');