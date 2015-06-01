function [cones_set,data_stuffing_set,dims_set]=codegenerate(LL, K, N1)
%Input:
% LL: # RAUs
% K: # MUs
% N1: # antennas at each RAU

%Output:
%cones_set: descriptions of the cone for each standard form
%data_stuffing_set: data (A, b, c) for each standard form
%dims_set: descriptions of the dimensions for each standard form

%%Network Setting%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

[cones_temp,data_stuffing_temp, dims_temp] = prob_to_socp_GSBF(params); % generate all the problem structures

cones_set=[cones_set,cones_temp];
data_stuffing_set=[data_stuffing_set, data_stuffing_temp];
dims_set=[dims_set,dims_temp];

end

end
