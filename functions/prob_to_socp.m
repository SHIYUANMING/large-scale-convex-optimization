function [data_socp,cones, dims, data_stuffing] = prob_to_socp(params)
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


%%%%%%%%%Matrix Stuffing+SCS%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Standard SOCP Data%%%%%%%%%%%%%%%%%%%%
B=params.K*(2*sum(params.N_set));  %total # beamforming coefficients ('2' is due to complex variables)

n=1+L+K+B;  
m=(L+K)+(B+1)+sum(K.*2*params.N_set+1)+K*(2*K+2); %dimensions of the cone

c=zeros(n,1); %objective vector
b=zeros(m,1);

cones.f=0;  %dimensions of cones
cones.l=L+K;   
cones.q=[B+1;K.*2*params.N_set+1; (2*K+2)*ones(K,1)]; 
cones.s=[];
 
Fi=[]; Fj=[]; Fv=[]; %store the index and value for the fixed part of the SOCP data
Qi=[]; Qj=[]; Qv=[]; %store the index and value for the problem data in SOCP data for matrix stuffing
Ci=[]; Cj=[]; Cv=[]; %store the index and value for the problem data in SOCP data for matrix stuffing

%%Stuffing the Objective Vector%%%%%%%%%%%%%
    c(1:1) = 1;

%%ADMM-Compliant Form for Objective Function%%
%%
Fi=[Fi;L+K+1]; Fj=[Fj;1]; Fv=[Fv;-1];
Fi=[Fi;(L+K+2:n)']; Fj=[Fj;(L+K+2:n)']; Fv=[Fv;-ones(B,1)];

%%ADMM-Compliant Form for Per-RRH Transmit Power Constraints%%
for l=1:L
%for linear constraints
Fi=[Fi;l]; Fj=[Fj;l+1]; Fv=[Fv;1]; b(l)=sqrt(params.P_set(l));  
   
%for SOC cone
D_temp=zeros(2*params.N_set(l), 2*sum(params.N_set)); 
if l==1
    Fi=[Fi;L+K+B+2]; Fj=[Fj;2]; Fv=[Fv;-1];
    
    D_temp([1:params.N_set(1)],[1:params.N_set(1)])=eye(params.N_set(1));
    D_temp([params.N_set(1)+1:2*params.N_set(1)], [sum(params.N_set)+1:sum(params.N_set)+sum(params.N_set(1))])=eye(params.N_set(1));
    
    D_temp=kron(eye(K), D_temp);
    D_temp=sparse(D_temp);
    
    Fi=[Fi;(L+K+B+2)+(1+mod(find(D_temp)-1,size(D_temp,1)))];
    Fj=[Fj;(L+K+1)+(1+floor((find(D_temp)-1)/size(D_temp,1)))];
    Fv=[Fv;nonzeros(-D_temp)];
  
else
    Fi=[Fi;L+K+(B+1)+sum(K*2*N_set(1:l-1)+1)+1]; Fj=[Fj;l+1]; Fv=[Fv;-1];
    
    D_temp([1:params.N_set(l)], [sum(N_set(1:l-1))+1:sum(N_set(1:l))])=eye(N_set(l));  
    D_temp([params.N_set(l)+1:2*params.N_set(l)], [sum(params.N_set)+sum(N_set(1:l-1))+1:sum(params.N_set)+sum(N_set(1:l))])=eye(N_set(l)); 
     
    D_temp=kron(eye(K), D_temp);
    D_temp=sparse(D_temp);
    
    Fi=[Fi;(L+K+(B+1)+sum(K*2*N_set(1:l-1)+1)+1)+(1+mod(find(D_temp)-1,size(D_temp,1)))];
    Fj=[Fj;(L+K+1)+(1+floor((find(D_temp)-1)/size(D_temp,1)))];
    Fv=[Fv;nonzeros(-D_temp)];
end
end

%%ADMM-Compliant Form for QoS Constraints%%%%%%%%%%
%%

for k=1:K
    %for linear constraints
    r_temp=zeros(B,1);
    if k==1
        r_temp(1:2*sum(params.N_set))=[real(params.H(:,k));imag(params.H(:,k))];
    else
        r_temp((k-1)*2*sum(params.N_set)+1: k*2*sum(params.N_set))=[real(params.H(:,k));imag(params.H(:,k))];
    end
   % r_temp=(sqrt(1+1/params.r_set(k)).*r_temp);
    r_temp=sparse(sqrt(1+1/params.r_set(k)).*r_temp);
   
    Fi=[Fi;L+k]; Fj=[Fj;L+k+1]; Fv=[Fv;1];  %fixed part
   
    Qi=[Qi;L+(k-1)+(1+floor((find(r_temp)-1)/size(r_temp,1)))];   %SOCP data dependents on SINR thresholds
    Qj=[Qj;(1+L+K)+(1+mod(find(r_temp)-1,size(r_temp,1)))]; 
    Qv=[Qv;nonzeros(-r_temp)];
    
    %for SOC cone constraints
    b(L+K+B+1+sum(K*2*params.N_set+1)+(k-1)*(2*K+2)+2*K+2)=params.delta_set(k);
    
    Fi=[Fi;(L+K+B+1)+sum(K*2*params.N_set+1)+(k-1)*(2*K+2)+1]; Fj=[Fj;1+L+k]; Fv=[Fv;-1];
    
    C_temp=kron(eye(K), [real(params.H(:,k)),-imag(params.H(:,k)); imag(params.H(:,k)), real(params.H(:,k))]'); 
    C_temp=sparse(C_temp); 
 
    Ci=[Ci;(L+K+B+1)+(sum(K*2*params.N_set+1)+(k-1)*(2*K+2)+1)+(1+mod(find(C_temp)-1,size(C_temp,1)))];
    Cj=[Cj;(L+K+1)+(1+floor((find(C_temp)-1)/size(C_temp,1)))];
    Cv=[Cv;nonzeros(-C_temp)];
   
end


%%Convert from sparse triplet to column compressed format.
%%
Ai=[Qi;Ci;Fi]; Aj=[Qj;Cj;Fj]; Av=[Qv;Cv;Fv];
A=sparse(Ai, Aj, Av, m, n);

%build data
data_socp=struct('c', c, 'b', b, 'A', A);
data_stuffing=struct('c', c, 'b', b,'Qv', Qv, 'Cv', Cv, 'Fv', Fv,'Ai', Ai, 'Aj', Aj, 'n',n,'m',m);
dims.n=n; dims.B=B;
