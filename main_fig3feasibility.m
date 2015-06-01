%This Matlab script can be used to generate the results by "Matrix Stuffing+SCS" in Fig.3 
%for infeasibility detection in the paper:

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
L=50; K=50; N1=1; Area=4*10^3; N_set=N1*ones(L,1); %'L': # RAU; 'K': # MUs; 'N1: # antennas in each RAU'
Q=[4:0.5:11]';  %QoS in dB
LC=1; % # loops for channel realizarions
prob=zeros(3,1); 
 %"prob(1): stuffing channel realization"
 %"prob(2):stuffing sinr thresholds"
 %"prob(3): stuffing transmit power"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
params.L=L;   %'L': # RRHs
params.K=K;    %'K': # MUs
params.N_set=N_set;  %set of antennas at all the RRHs
params.delta_set=1*ones(K,1); %set of noise covariance
params.P_set=10^(0)*ones(L,1);   %set of transmit power constraints at all the RRHs


%%%%%%%%%Generate SOCP Structure%%%%%%%%%%%%%
%%
%can be done offline
params.P_set=ones(params.L,1);
params.r_set=ones(params.K,1);  
params.H=ones(sum(params.N_set), params.K)+i*ones(sum(params.N_set), params.K);
[data_socp,cones,dims,data_stuffing] = prob_to_socp(params); % generate and store standard form

%%%%%%%%%%Parameters for SCS solver%%%%%%%%%%%%
%%
paramscs.VERBOSE = 0;paramscs.EPS = 1.00e-03; paramscs.ALPHA = 1.8;
paramscs.MAX_ITERS = 10^4; 
paramscs.NORMALIZE = 1; paramscs.SCALE =20;

%%%%%%initialize results%%%%%%%%%%%%%%%%%%%%%%
%%
SU_matscs=0;  SU_matscs_temp=0;%count the success of channel realizations
TransTime_temp=0; % modeling time
SolveTime_temp=0; % solving time
TotalTime_temp=0; % total time

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
tic;
prob=[1;1;1];
[data_socp, data_stuffing] = matrix_stuffing(params,data_stuffing,prob); %standard form data stuffing
TransTime(lq)=toc;  % modeling time

%call SCS solver
tic;
if lq>1
    %%Warm-starting
    data_socp.x = xx; data_socp.y = yy; data_socp.s = ss;
end
[xx,yy,ss, info] = scs_direct(data_socp,cones,paramscs); %call SCS solver

SolveTime(lq)=toc; % solving time

if strfind(info.status,'Solved') 
    SU_matscs(lq)=1; % feasible 
else
    SU_matscs(lq)=0; % infeasible
end

    end
 
SU_matscs_temp=SU_matscs_temp+SU_matscs; 
TransTime_temp=TransTime_temp+TransTime;
SolveTime_temp=SolveTime_temp+SolveTime;
TotalTime_temp=TotalTime_temp+(TransTime+SolveTime);
    
end

%%Prepare outputs
SU_matscs=SU_matscs_temp/LC; % feasiblity 
TransTime=TransTime_temp/LC; % modeling time
SolveTime=SolveTime_temp/LC; % solving time
TotalTime=TotalTime_temp/LC; % total time
    

plot(Q, SU_matscs*100, 'b-d','LineWidth',2.0, 'MarkerSize',8); %CVX
hold on;
h=legend('Matrix Stuffing+SCS');

%
xlabel('Target SINR [dB]','fontsize',14,'fontweight','b','fontname','helvetica');
ylabel('Empirical probability of feasibility','fontsize',14,'fontweight','b','fontname','helvetica');

h=gca;
labels=get(h,'yticklabel'); 
for i=1:size(labels,1)
   labels_modif(i,:)=[labels(i,:) '%']; 
end
set(h,'yticklabel',labels_modif);
