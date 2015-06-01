function [H] =channel_realization(L, K, N_set, Area)
% Generate channel matrix

%INPUT:
%L      =# RRHs
%K      =# MUs
%N_set    =# antennas at each RRH
%Area      =length of the area

%OUTPUT:
% H     = (sum(N_set) x K) channel matrix

%%%%%%%%%%%%%%%%Network Realization%%%%%%%%%%%%%%%%%%%%%%%%
 U_position=Area.*(rand(2,K)-0.5);  %% user positions
 B_position=Area.*(rand(2,L)-0.5);  %%RRH positions

 %%%%%Generate Large-Scale Fading%%%%%%%%%%%%%
for k=1:K
    for l=1:L
                   d=(norm(B_position(:,l)-U_position(:,k))+10);
                   D(l,k)=4.4*10^(5)/(d^(1.88)*10^(normrnd(0,6.3)/20));
    end
end

%%%%%%Generate Small-Scale Fading%%%%%%%%%%%%%
for k=1:K
    for l=1:L
        if l==1
            H(1:N_set(1),k)=D(l,k)*(normrnd(0,1/sqrt(2),N_set(1),1)+i*normrnd(0,1/sqrt(2),N_set(1),1));  %%%nosie normalized to 1
        else
           H(sum(N_set(1:(l-1)))+1:sum(N_set(1:l)),k)=D(l,k)*(normrnd(0,1/sqrt(2),N_set(l),1)+i*normrnd(0,1/sqrt(2),N_set(l),1));  %%%nosie normalized to 1
        end
    end
end
