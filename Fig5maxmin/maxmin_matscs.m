function [r, TransTime, SolveTime] = maxmin_matscs(params,data_stuffing, cones)

%%%%%%%%%%Parameters for SCS Solver%%%%%%%%%%%%
paramscs.VERBOSE = 0;
paramscs.EPS = 1.00e-03; paramscs.ALPHA = 1.80;
paramscs.MAX_ITERS = 10^4; paramscs.NORMALIZE = 1; 
paramscs.SCALE = 20;

%%%%%%%%Maxmin Algorithm%%%%%%%%%%%%%%
K=params.K;
r_low=params.r_low; r_up=params.r_up;  epsilon=params.epsilon; 

prob=[1;1;1]; %stuff for SINR thresholds

LP_counter=0;

TransTime=0; SolveTime=0;
tic;

while r_up-r_low>epsilon
    LP_counter=LP_counter+1;
    
r=(r_low+r_up)/2;
params.r_set=r*ones(params.K,1);  %set of SINR thresholds

%% socp data stuffing   
tic;
[data_socp, data_stuffing] = matrix_stuffing(params,data_stuffing,prob);
TransTime=TransTime+toc;

%% call SCS solver
if LP_counter>1
    %%Warm-starting
    data_socp.x = xx; data_socp.y = yy; data_socp.s = zz;
end

tic;
[xx,yy, zz, info] = scs_direct(data_socp,cones,paramscs);
SolveTime=SolveTime+toc;


if strfind(info.status,'Solved') 
    r_low=r;
else
    r_up=r;
end
end

r=(r_up+r_low)/2;