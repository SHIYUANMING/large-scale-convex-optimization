clc;clear all;
%cvx_solver sedumi
%cvx_solver sdpt3
%cvx_solver mosek

cvx_solver scs
cvx_solver_settings('MAX_ITERS', 10^4, 'SCALE', 20);
%cvx_quiet(true)

%addpath('/Users/yuanmingshi/scs/');
%addpath('/Users/yuanmingshi/Dropbox/Research/Software/pdos-master/matlab/');
%cvx_solver pdos

%%Problem Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=100; K=100; N1=1; Area=4*10^3; N_set=N1*ones(L,1);
Q=[4:0.5:11]';  %QoS in dB
LC=200; % # loops for channel realizarions
savedata=true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.L=L;   %'L': # RRHs
params.K=K;    %'K': # MUs
params.N_set=N_set;  %set of antennas at all the RRHs
params.delta_set=1*ones(K,1); %set of noise covariance
params.P_set=10^(0)*ones(L,1);   %set of transmit power constraints at all the RRHs

%%%%%%initialize results%%%%%%%%%%%%%%%%%%%%%%
SU_CVX=0;  SU_CVX_temp=0;   %count the success of channel realizations
SolveTime_CVX=0;  SolveTime_CVX_temp=0;  %Time for solver in CVX
TotalTime_CVX=0;  TotalTime_CVX_temp=0;  %Time for solver in CVX

for ss=1:LC   %generate channel
H(:,:,ss)=channel_realization(L, K, N_set, Area);
end

% load('H.mat');

for lp=1:LC
    params.H=H(:,:,lp);
    
    for lq=1:length(Q)
        params.r_set=10^(Q(lq)/10)*ones(K,1);  %set of SINR thresholds   
     
%%%%%%%%%%CVX%%%%%%%%%%%%%%%%%%%%%%%
%Use Software CVX
tic;
[Wsolution_cvx, feasible_cvx, solving_time] = powermin_cvx(params);
TotalTime_CVX(lq)=toc;

%%%%%%%%%%%%%%%%%%%%
SU_CVX(lq)=feasible_cvx;
SolveTime_CVX(lq)=solving_time;
    
    end 
    SU_CVX_temp=SU_CVX_temp+SU_CVX;
    SolveTime_CVX_temp=SolveTime_CVX_temp+SolveTime_CVX;
    TotalTime_CVX_temp=TotalTime_CVX_temp+TotalTime_CVX;
end

%% Final Results
TotalTime_CVX=TotalTime_CVX_temp/LC;
SolveTime_CVX=SolveTime_CVX_temp/LC;
SU_CVX=SU_CVX_temp/LC;

if savedata==1
save('TotalTime_CVX.mat','TotalTime_CVX');
save('SolveTime_CVX.mat','SolveTime_CVX');
save('SU_CVX.mat','SU_CVX');
end

plot(Q, SU_CVX*100, 'g-s','LineWidth',2.0, 'MarkerSize',8); %CVX
hold on;

h=legend('CVX');
xlabel('Target SINR [dB]','fontsize',14,'fontweight','b','fontname','helvetica');
ylabel('Empirical probability of feasibility','fontsize',14,'fontweight','b','fontname','helvetica');

h=gca;
labels=get(h,'yticklabel'); 
for i=1:size(labels,1)
   labels_modif(i,:)=[labels(i,:) '%']; 
end
set(h,'yticklabel',labels_modif);

%MySendMail;