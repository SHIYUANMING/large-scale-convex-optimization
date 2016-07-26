function [data_socp] = matrix_stuffing(params,data_stuffing)
% matrix stuffing: copy problem data to the socp data
% params: problem data
%data_stuffing: socp data structure from function prob_to_socp
%stuffing.b: vector b
%stuffing.Q: depends on SINR thresholds and channel realizations
%stuffing.C: depends on channel realizations
%prob (3x1) vector: "prob(1): stuffing channel realization", "prob(2):
%stuffing sinr thresholds", "prob(3): stuffing transmit power".

prob=[1;1;1];

H=params.H;
%%%%%%%Matrix Stuffing%%%%%%%%%%%%%%%%%%%
%for stuffing.Q
if prob(1)==1|prob(2)==1
H=[real(params.H);imag(H)];
data_stuffing.Qv=-vec(repmat(sqrt(1+(1./params.r_set)'), sum(2*params.N_set),1).*H);
end

%for stuffing.C
if prob(1)==1
    H_temp=zeros(4*sum(params.N_set),params.K);
    H_temp(1:2:2*sum(params.N_set),:)=real(params.H);
    H_temp(2:2:2*sum(params.N_set),:)=-imag(params.H);
    H_temp(2*sum(params.N_set)+1:2:4*sum(params.N_set),:)=imag(params.H);
    H_temp(2*sum(params.N_set)+2:2:4*sum(params.N_set),:)=real(params.H);
    
    data_stuffing.Cv=-vec(repmat(H_temp,params.K,1));
end

%for stuffing.b
if prob(3)==1
data_stuffing.b([1:params.L])=sqrt(params.P_set);
end

%%Convert from sparse triplet to column compressed format.
%%
Av=[data_stuffing.Qv;data_stuffing.Cv;data_stuffing.Fv];
A=sparse(data_stuffing.Ai, data_stuffing.Aj, Av, data_stuffing.m, data_stuffing.n);

%build data
data_socp=struct('c', data_stuffing.c, 'b', data_stuffing.b, 'A', A);
data_stuffing_new=struct('c', data_stuffing.c, 'b', data_stuffing.b,'Qv', data_stuffing.Qv, 'Cv', data_stuffing.Cv, 'Fv', data_stuffing.Fv,'Ai', data_stuffing.Ai, 'Aj', data_stuffing.Aj, 'n',data_stuffing.n,'m',data_stuffing.m);
