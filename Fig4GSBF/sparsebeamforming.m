function [feasible,z_sparse, Wsolution] = sparsebeamforming(params)

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
weights=params.weights;

%%%%%%%%CVX+SCS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time_pat_cvx = 'Total CPU time \(secs\)\s*=\s*(?<total>[\d\.]+)';

%tic
cvx_begin quiet
variable W(sum(params.N_set), K) complex;   %Variable for N x K beamforming matrix
variable t(L,1);   % Auxiliary variable for power constraints
minimize (weights*t)
subject to
     for l=1:L    %%%Transmit Power Constraints
         %%
         if l==1
             norm(W(1:params.N_set(1),:),'fro')<=t(l);
             t(l)<=sqrt(params.P_set(l));
         else
             norm(W(sum(params.N_set(1:l-1))+1:sum(params.N_set(1:l)),:),'fro')<=t(l);
             t(l)<=sqrt(params.P_set(l));
         end
     end
     
     for k=1:K        %%%%%%%%%QoS constraints
      %%
         norm([params.H(:,k)'*W, params.delta_set(k)],'fro')<=sqrt(1+1/params.r_set(k))*real(params.H(:,k)'*W(:,k));   
     end
     cvx_end
%output = evalc('cvx_end');
 
%build output
%%
     if  strfind(cvx_status,'Solved') 
         feasible=true;
         z_sparse=t;
         Wsolution=W;
     else
         feasible=false;
         z_sparse=[];
         Wsolution=[];
     end