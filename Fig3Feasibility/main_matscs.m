clc;clear all;
cvx_solver scs

 addpath('/home/yuanming/Documents/Software/scs-master/matlab');
%addpath('/Users/yuanmingshi/Dropbox/Research/Software/scs-master/matlab');

%%Problem Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=100; K=100; N1=1; Area=4*10^3; N_set=N1*ones(L,1);
%Q=2;
Q=[4:0.5:11]';  %QoS in dB
%Q=10;
LC=200; % # loops for channel realizarions
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
[data_socp,cones,dims,data_stuffing] = prob_to_socp(params);

%%%%%%%%%%Parameters for SCS solver%%%%%%%%%%%%
%%
paramscs.VERBOSE = 0;paramscs.EPS = 1.00e-03; paramscs.ALPHA = 1.8;
paramscs.MAX_ITERS = 10^4; 
paramscs.NORMALIZE = 1; paramscs.SCALE =20;

%%%%%%initialize results%%%%%%%%%%%%%%%%%%%%%%
%%
SU_matscs=0;  SU_matscs_temp=0;%count the success of channel realizations
TransTime_temp=0;
SolveTime_temp=0;
TotalTime_temp=0;

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
tic;
prob=[1;1;1];
[data_socp, data_stuffing] = matrix_stuffing(params,data_stuffing,prob);
TransTime(lq)=toc;

%call SCS solver
%%
tic;
if lq>1
    %%Warm-starting
    data_socp.x = xx; data_socp.y = yy; data_socp.s = ss;
end
[xx,yy,ss, info] = scs_direct(data_socp,cones,paramscs);

SolveTime(lq)=toc;

if strfind(info.status,'Solved') 
    SU_matscs(lq)=1;
else
    SU_matscs(lq)=0;
end

    end
 
SU_matscs_temp=SU_matscs_temp+SU_matscs;
TransTime_temp=TransTime_temp+TransTime;
SolveTime_temp=SolveTime_temp+SolveTime;
TotalTime_temp=TotalTime_temp+(TransTime+SolveTime);
    
end

%%Final Results
SU_matscs=SU_matscs_temp/LC;
TransTime=TransTime_temp/LC;
SolveTime=SolveTime_temp/LC;
TotalTime=TotalTime_temp/LC;

    save('SU_matscs.mat','SU_matscs');
    save('TransTime.mat','TransTime');
    save('SolveTime.mat','SolveTime');
    save('TotalTime.mat','TotalTime');
    


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

%MySendMail;
