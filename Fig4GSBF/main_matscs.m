clc;clear all;
%cvx_solver sedumi
%cvx_solver sdpt3
%cvx_solver mosek

cvx_solver scs
cvx_solver_settings('SCALE',20, 'MAX_ITERS', 10^4);

addpath('/Users/yuanmingshi/Dropbox/Research/Software/scs-master/matlab');

%%Problem Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=50; K=50; N1=2; Area=6*10^3; N_set=N1*ones(L,1);

Pc=5.6+0*[1:1:L];    %fronthaul network power consumption
coeff=4*ones(1,L);  %power amplifer coefficients

Q=[0:2:6]';  %QoS in dB
%Q=4;
LC=50; %# loops for channel realizarions

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

%load('H.mat');

load('cones_set.mat');
load('data_stuffing_set.mat');
load('dims_set');

for lp=1:LC
    H_orignal=H(:,:,lp);
    
    for lq=1:length(Q)
        
        TransTime_record=0;  %record time at each SINR target
        SolveTime_record=0;
       
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
         [data_socp] = matrix_stuffing(params,data_stuffing_set(params.L));
         TransTime_record=TransTime_record+toc;
         
         tic;
         [xx,yy,ss, info] = scs_direct(data_socp,cones_set(params.L),paramscs);
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
A_number=A_number_temp./SU_counter;
TotalPower=TotalPower_temp./SU_counter;
TransmitPower=TransmitPower_temp./SU_counter;

TransTime=TransTime_temp/SU_counter;
SolveTime=SolveTime_temp/SU_counter;
TotalTime=TotalTime_temp/SU_counter;

save('A_number.mat','A_number');
save('TotalPower.mat','TotalPower');
save('TransmitPower.mat','TransmitPower');

save('TransTime.mat','TransTime');
save('SolveTime.mat','SolveTime');
save('TotalTime.mat','TotalTime');

plot(Q,TotalPower,'b-d','LineWidth',1.2, 'MarkerSize',8); %Proposed Bi-section GSBF
hold on;

h=legend('GSBF');
xlabel('Target SINR [dB]','fontsize',12,'fontweight','b','fontname','helvetica');
ylabel('Average Network Power Consumption [W]','fontsize',14,'fontweight','b','fontname','helvetica');