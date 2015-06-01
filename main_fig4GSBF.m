%This Matlab script can be used to generate the results by "Matrix Stuffing+SCS" in Fig.4 in the
%paper:

%Yuanming Shi, Jun Zhang, Brendan O'Donoghue, and Khaled B. Letaief, "Large-scale convex 
%optimization for dense wireless cooperative networks," IEEE Trans. Signal Process., to appear, 2015.

% In particular, Fig.4 is an example the "bi-section group sparse beamforming (GSBF)
% algorithm" in the paper:

%Yuanming Shi, Jun Zhang, and Khaled B. Letaief, "Group sparse beamforming for green 
%Cloud-RAN," IEEE Trans. Wireless Commun., vol. 13, pp. 2809? 2823, May 2014.

%The implementation of this algorithm needs to download the solver "scs" 
% https://github.com/cvxgrp/scs 
%based on the paper

%Brendan O'Donoghue, E. Chu, N. Parikh, and S. Boyd, "Conic optimization via operator splitting and 
%homogeneous self-dual embedding," arXiv preprint arXiv:1312.3039, 2013.

%This is version 1.0 (Last edited: 2015-06-01)

%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.

%Note: the channels are generated randomly, thus the results
%will not be exactly the same as in the paper.

clc;clear all;

%add the path for the solver "scs"
%addpath('/scs-master/matlab');
addpath('/Users/Yuanming/Dropbox/Research/Software/SCS/scs-master/matlab');

%%Problem Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=20; K=20; N1=2; Area=2*10^3; N_set=N1*ones(L,1); %'L': # RAU; 'K': # MUs; 'N1: # antennas in each RAU'

Pc=5.6+0*[1:1:L];    %fronthaul network power consumption
coeff=4*ones(1,L);  %power amplifer coefficients

Q=[0:2:6]';  %QoS in dB
LC=5; %# loops for channel realizarions

prob=[1;1;1]; %"prob(1): stuffing channel realization"  %"prob(2):stuffing sinr thresholds" %"prob(3): stuffing transmit power"
%%%%%%%%%%Parameters for SCS Solver%%%%%%%%%%%%
paramscs.VERBOSE = 0;paramscs.EPS = 1.00e-03; paramscs.ALPHA = 1.8;
paramscs.MAX_ITERS = 10^4; 
paramscs.NORMALIZE = 1; paramscs.SCALE = 20;

%%%%%%initialize results%%%%%%%%%%%%%%%%%%%%%%
%%
SU_counter=0; %count the success of channel realizations
A_number_temp=0;  TotalPower_temp=0;  TransmitPower_temp=0;
TransTime_temp=0; SolveTime_temp=0; TotalTime_temp=0;
% 
for ss=1:LC   %generate channel
H(:,:,ss)=channel_realization(L, K, N_set, Area);
end

%%Generate and store the standard form problem structures for all the
%%candidate network sizes: this can be done offline

[cones_set,data_stuffing_set,dims_set]=codegenerate(L, K, N1);

%cones_set: descriptions of the cone for each standard form
%data_stuffing_set: data (A, b, c) for each standard form
%dims_set: descriptions of the dimensions for each standard form

for lp=1:LC
    H_orignal=H(:,:,lp);
    
    for lq=1:length(Q)
        
        TransTime_record=0;  %record time at each SINR target for modeling time
        SolveTime_record=0;   %record time at each SINR target for solving time
       
        %%%%%%%%%%%Problem Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        params.H=H_orignal;
        params.L=L;   %'L': # RRHs
        params.N_set=N_set;  %set of antennas at all the RRHs
        params.P_set=10^(0)*ones(L,1);   %set of transmit power constraints at all the RRHs
        
        params.K=K;    %'K': # MUs
        params.delta_set=1*ones(K,1); %set of noise covariance
        params.r_set=10^(Q(lq)/10)*ones(K,1);  %set of SINR thresholds   
        
%%%%%%%%%%Group Sparse Beamforming%%%%%%%%%%%%%%%%%%%%%%%
%%L_l2 norm minimization
%%
weights = sqrt(Pc);  P_optimal=10^100; 

[feasible_cvx,z_sparse] = sparsebeamforming(params, weights); %sparse beamforming solutions

%%
 if feasible_cvx==1 
     
     %%%RRH Ordering
     %%
     for l=1:L    
         if l==1
             g(l) = norm(params.H(1:N_set(1),:),'fro');
         else
             g(l) = norm(params.H(sum(N_set(1:l-1))+1:sum(N_set(1:l)),:),'fro');
         end
     end
     
     SparseW=sqrt(Pc)./g; %1 x L Group sparse ording weights
     
     [BS,BS_index]=sort(z_sparse'./SparseW);   %RRH ordering
     
     D_set=[]; A_set=BS_index;  %A_set: active RRH set, D_set: inactive RRH set
     D_optimal=[]; A_optimal=BS_index;  % optimal active and inactive RRH set 
     
     %%%%%%%%%%%%%%%%%%Process Deflation Procedure%%%%%%%%%%%%%%%%%%%%%%
     %%
     for A=0:L-2
         
         %% Updata Problem Data 
         params.L=length(A_set);   %'L': # RRHs
         params.N_set=N1*ones(params.L,1);  %set of antennas at all the RRHs
         params.P_set=10^(0)*ones(params.L,1);   %set of transmit power constraints at all the RRHs
         
         %% Updata Channel
         H_index=N1*A_set;
         if N1>1
         for nn=1:(N1-1)
             H_index=[H_index, N1*A_set-nn*ones(length(A_set),1)'];
         end
         end
         params.H = H_orignal(H_index,:);
         
         
         %% Call SCS with Matrix Stuffing
         tic;
         [data_socp] = matrix_stuffing(params,data_stuffing_set(params.L),prob); %standard form data stuffing
         TransTime_record=TransTime_record+toc; 
         
         tic;
         [xx,yy,ss, info] = scs_direct(data_socp,cones_set(params.L),paramscs); % call SCS solver
         SolveTime_record=SolveTime_record+toc;
         
         feasible=strfind(info.status,'Solved');
         
         if feasible==1
             P_optimal=coeff(1)*norm(xx(dims_set(params.L).n-dims_set(params.L).B+1:dims_set(params.L).n),'fro')^2+sum(Pc(A_set));  %network power consumption
             D_optimal=D_set; A_optimal=A_set;
             D_set=[D_set, BS_index(A+1)];  %%%process next selection: updata active RRH set and number
             A_set=setdiff([1:L],D_set); 
         elseif feasible==0&isempty(D_set)==1   %infeasible
             break;  
         else feasible==0&isempty(D_set)==0   %infeasible
             D_set=setdiff(D_set,D_set(length(D_set)));   %%%%%%Backup one step to recoder the final feasible set
             A_set=setdiff([1:L], D_set);
             break,
         end
     end
 end
 
 %% Results for each SINR
 TotalPower(lq)=P_optimal;
 A_number(lq)=size(A_optimal,2);
 TransmitPower(lq)=P_optimal-sum(Pc(A_optimal));
 
 TransTime(lq)=TransTime_record;
 SolveTime(lq)=SolveTime_record;
 
    
   end
 
%%
    if norm(TotalPower)<10^100
     SU_counter=SU_counter+1;
     A_number_temp=A_number_temp+A_number;
     TotalPower_temp=TotalPower_temp+TotalPower;
     TransmitPower_temp=TransmitPower_temp+TransmitPower;
     
     TransTime_temp=TransTime_temp+TransTime;
     SolveTime_temp=SolveTime_temp+SolveTime;
     TotalTime_temp=TotalTime_temp+(TransTime+SolveTime);
     
    end
end

%% Final average results
A_number=A_number_temp./SU_counter;  % active RAUs
TotalPower=TotalPower_temp./SU_counter; % total power consumption
TransmitPower=TransmitPower_temp./SU_counter; % transmit power consumption

TransTime=TransTime_temp/SU_counter; % modeling time
SolveTime=SolveTime_temp/SU_counter; % solving time
TotalTime=TotalTime_temp/SU_counter; % total time

plot(Q,TotalPower,'b-d','LineWidth',1.2, 'MarkerSize',8); %Proposed Bi-section GSBF
hold on;

h=legend('GSBF');
xlabel('Target SINR [dB]','fontsize',12,'fontweight','b','fontname','helvetica');
ylabel('Average Network Power Consumption [W]','fontsize',14,'fontweight','b','fontname','helvetica');