clc;clear all;
cvx_solver scs

addpath('/Users/yuanmingshi/Dropbox/Research/Software/scs-master/matlab');

%%Network Setting%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL=50; K=50; N1=2;
L_set=[1:1:LL]; 

cones_set=[];
data_stuffing_set=[];
dims_set=[];

for l=1:length(L_set)
%% Problem Data
L=L_set(l);
N_set=N1*ones(L,1);

params.L=L;   %'L': # RRHs
params.K=K;    %'K': # MUs
params.N_set=N_set;  %set of antennas at all the RRHs
params.delta_set=1*ones(K,1); %set of noise covariance
params.P_set=10^(0)*ones(L,1);   %set of transmit power constraints at all the RRHs


%% Generate SOCP Structure
params.P_set=ones(params.L,1);
params.r_set=ones(params.K,1);  
params.H=ones(sum(params.N_set), params.K)+i*ones(sum(params.N_set), params.K);

[cones_temp,data_stuffing_temp, dims_temp] = prob_to_socp(params);

cones_set=[cones_set,cones_temp];
data_stuffing_set=[data_stuffing_set, data_stuffing_temp];
dims_set=[dims_set,dims_temp];

end

save('cones_set.mat', 'cones_set');
save('data_stuffing_set.mat', 'data_stuffing_set', '-v7.3');
save('dims_set.mat', 'dims_set');

% %socp data stuffing   
% %%
% prob=[1;1;1];
% 
% [data_socp] = matrix_stuffing(params,data_stuffing,prob);
% 
% [xx,yy,ss, info] = scs_direct(data_socp,cones,paramscs);
