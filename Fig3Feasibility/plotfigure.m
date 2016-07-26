load('SU_matscs.mat');
load('SU_scs.mat');
load('SU_mosek.mat');


Q=[2:0.5:8]';  %QoS in dB
plot(Q, SU_matscs*100, 'gd','LineWidth',2.0, 'MarkerSize',12); %CVX
hold on;
plot(Q, SU_scs*100, 'ro','LineWidth',2.0, 'MarkerSize',8); %CVX
hold on;
plot(Q, SU_mosek*100, 'b-.','LineWidth',2.0, 'MarkerSize',8); %CVX
hold on;
h=legend('Matrix Stuffing+SCS', 'CVX+SCS', 'CVX+MOSEK');

xlabel('Target SINR [dB]','fontsize',14,'fontweight','b','fontname','helvetica');
ylabel('Empirical probability of feasibility','fontsize',14,'fontweight','b','fontname','helvetica');