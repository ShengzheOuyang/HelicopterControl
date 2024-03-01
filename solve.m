syms F_f F_b
x=sym("x",[6,1]);
dx = sym("dx",[6,1]);

%% init parameters
% parameters definetion
l1 = 0.655;            %length of main arm
l2 = 0.262;            % length of Secondary arm [m]
l3 = 0.1355;           % length of midpoint of Secondary arm to crossing point[m]
l4 = 0.2285;           % length of joint 1 to crossing point
l5 = 0.21;             % length of Joint 1 to midpoint of main arm [m]
l6= 0.1775;            % half of the distance between the motors
m = 1.87;             % weight of counterweight [kg]
m_arms = 0.138;        % weight of Secondary arm [kg]
m_arml = 0.377;        % weight of Main arm [kg]
mb = 0.661;        % weight of Back motor + Bodies [kg]
mf = 0.661;        % weight of Front motor + Bodies [kg]
g = 9.81;               % gravity [m/s^2]
% moment of inertia for beta
J_beta = (mf+mb)*(l1^2) ...,
        +m*(l2+l4)^2 +m_arms*(1/12*l2^2+(l3+l4))^2 ...,
        +m_arml*(1/12*(l1+l4)^2+l5^2);

% moment of inertia for alpha
J_alpha = (mf + mb)*l1^2 + m*(l2 + l4)^2 ...,
    + m_arms*(1/12*l2^2+(l3+l4)^2) ...,
    + m_arml*(1/12*(l1+l4)^2+l5^2) + ...,
    (mf+mb) * l6^2;

% moment of inertia for gamma
J_gamma = (mf+mb)*(l6^2);

% state space
dx(1) = x(4);
dx(2) = x(5);
dx(3) = x(6);
dx(4) = - (F_f + F_b)*sin(x(3))*l1/J_alpha;
dx(5) = ((F_f + F_b)*l1+ m*g*(l2+l4)-m_arms*g*(l3+l4)-m_arml*g*l5-(mf+mb)*g*l1 ) / J_beta;
dx(6) = ((F_f - F_b)*l6)/J_gamma;

% subs(dx,x,[0;0;-5;0;0;0]);
A_= jacobian(dx,x);
B_= jacobian(dx,[F_f,F_b]);
syms x1 x2 x3 x4 x5 x6
A=eval(subs(A_,{x1,x2,x3,x4,x5,x6,F_f,F_b},{0,-5*pi/180,0,0,0,0, u_ap(1)   u_ap(2)}));
B=eval(subs(B_,{x1,x2,x3,x4,x5,x6,F_f,F_b},{0,-5*pi/180,0,0,0,0, u_ap(1)   u_ap(2)}));

C=[1 0 0 0 0 0;
    0 1 0 0 0 0
    0 0 1 0 0 0];
    
% check Stability
pols = eig(A);
% check controllability and observability
beobachtbarkeit = rank(obsv(A,C)) == 6;
steubarkeit =  rank(ctrb(A,B)) == 6;
