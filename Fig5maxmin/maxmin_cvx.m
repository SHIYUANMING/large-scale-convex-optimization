function [r, TotalTime, SolveTime] = maxmin_cvx(params)

%prob_to_socp: maps PARAMS into a struct of SOCP matrices
%input struct 'parms' has the following fields:
%params.L;   %'L': # RRHs
%params.K;    %'K': # MUs
%params.N_set;  %set of antennas at all the RRHs
%params.delta_set; %set of noise covariance

%%%%%%%%%%%%%%Problem Instances%%%%%%%%%%%%%
%params.r_set;  %set of SINR thresholds
%params.H;  %Channel Realization
%params.P_set;   %set of transmit power constraints at all the RRHs

%%%%%%%%Problem Data%%%%%%%
K=params. K;
L=params. L;
N_set=params.N_set;


%%%%%%%%Maxmin Algorithm%%%%%%%%%%%%%%%%
r_low=params.r_low; r_up=params.r_up;  epsilon=params.epsilon; 

TotalTime=0; SolveTime=0;
tic;

while r_up-r_low>epsilon
r=(r_low+r_up)/2;
params.r_set=r*ones(params.K,1);  %set of SINR thresholds
%%%%%%%%CVX+SCS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cvx_begin 
variable W(sum(params.N_set), K) complex;   %Variable for N x K beamforming matrix
minimize norm(W,'fro')
subject to
     for l=1:L    %%%Transmit Power Constraints
         %%
         if l==1
             norm(W(1:params.N_set(1),:),'fro')<=sqrt(params.P_set(l));
         else
             norm(W(sum(params.N_set(1:l-1))+1:sum(params.N_set(1:l)),:),'fro')<=sqrt(params.P_set(l));
         end
     end

     for k=1:K        %%%%%%%%%QoS constraints
      %%
         norm([params.H(:,k)'*W, params.delta_set(k)])<=sqrt(1+1/params.r_set(k))*real(params.H(:,k)'*W(:,k));   
     end
%cvx_end

output = evalc('cvx_end');

%% Timing for SDPT3
%  time_pat_cvx = 'Total CPU time \(secs\)\s*=\s*(?<total>[\d\.]+)';  %SDPT3
%  timing = regexp(output, time_pat_cvx, 'names'); %SDPT3
%  solving_time=str2num(timing.total);

%% Timing for SCS
time_pat_cvx1 = 'Timing: Total solve time: (?<total>[\d\.e+\d]+)';  %SCS
time_pat_cvx2 = 'Timing: Total solve time: (?<total>[\d\.e-\d]+)';

timing = regexp(output, time_pat_cvx1, 'names');  %SCS
solving_time=str2num(timing.total);

if length(solving_time)==0
timing = regexp(output, time_pat_cvx2, 'names');
solving_time=str2num(timing.total);
end

%bi-section
%%
     if  strfind(cvx_status,'Solved') 
         r_low=r;
     else
         r_up=r;
     end

%% Timing 
SolveTime=SolveTime+solving_time;
end
TotalTime=toc;
r=(r_up+r_low)/2;