p_beo=[-25 -26 -27 -28 -29 -30];
L=place(A',C',p_beo)';
R=eye(2);
A_neu = [A,zeros(6,2);-C(1:2,:),zeros(2,2)];
B_neu = [B;zeros(2,2)];
C_neu = [C,zeros(3,2)];
q_neu = [1000; 500; 20; 10; 10; 1;20;20];
Q_neu = diag(q_neu);
[K_lqi,P_neu,e_neu] = lqr(A_neu,B_neu,Q_neu,R);
K = K_lqi(:,1:6);
V = -pinv(C*inv(A-B*K)*B);