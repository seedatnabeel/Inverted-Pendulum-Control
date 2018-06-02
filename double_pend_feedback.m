
%%%%%%%%%%%%%%%%%%%

% Define params for double pendulum
m0 = 1; 
m1 = 0.5; 
m2 = 0.25;
L1 = 0.25; 
L2 = 0.5;

% define for the G matrix
f1 = (m1/2 + m2)*L1*9.81; 
f2 = m2*L2*9.81/2;

%%%% matrices to linearize
D_=[m0+m1+m2 (m1/2+m2)*L1 (m2*L2)/2;(m1/2+m2)*L1 (m1/3+m2)*(L1*L1) (m2*L1*L2)/2; (m2*L2)/2 (m2*L1*L2)/2 (m2*L2*L2)/3];
G=[0 0 0; 0 -f1 0 ; 0 0 -f2];
H=[1; 0; 0 ];

% state space matrices for linear double inverted pendulum
A=[zeros(3,3) eye(3); -inv(D_)*G zeros(3,3)];
B=[zeros(3,1); inv(D_)*H];
C=[1 1 1 1 1 1];
D=0;

% find controllability matrix
controllability_matrix=horzcat(B, A*B, A^2*B, A^3*B, A^4*B, A^5*B);
double_ctrb_matrix=horzcat(B, A*B, A^2*B, A^3*B, A^4*B, A^5*B);

%check using the ctrb function
ctrb_func= ctrb(A,B);

unctrb_states=length(A) - rank(controllability_matrix); %test controllability
poles=eigs(A); % find the poles of A

% find observability matrix
double_obs_matrix=obsv(A,C);

% pole placement using matrix algebra
syms s k1 k2 k3 k4 k5 k6;
sI=[s 0 0 0 0 0 ; 0 s 0 0 0 0 ; 0 0 s 0 0 0; 0 0 0 s 0 0 ; 0 0 0 0 s 0 ; 0 0 0 0 0 s];
abs_poles=abs(poles);
desired_roots_non_repeated=roots([1,2*0.59115*3.383,3.383*3.383]); % non repeated roots
u1=desired_roots_non_repeated(1);
u2=desired_roots_non_repeated(2);
u3=desired_roots_non_repeated(1)*3;
u4=desired_roots_non_repeated(2)*3;
u5=desired_roots_non_repeated(1)*5;
u6=desired_roots_non_repeated(2)*5;
ctlr_poles=[u1 u2 u3 u4 u5 u6];
digits(5);
cl_char_eq=expand((s-u1)*(s-u2)*(s-u3)*(s-u4)*(s-u5)*(s-u6)); %closed loop characteristic equation

simplified_cl_char_eq=vpa(cl_char_eq);
K=[k1 k2 k3 k4 k5 k6];
lhs=vpa(det(sI-A+B*K)); % calculate the determinant 

pretty(simplified_cl_char_eq);
pretty(collect(lhs));

coeffs_rhs=coeffs(simplified_cl_char_eq);

eqn1 = 2217.3 * k1 == coeffs_rhs(1);
eqn2 = + 2217.3 * k4== coeffs_rhs(2);
eqn3 = (226.02 *k2 - 122.43* k1 + 226.02* k3 + 3880.2)== coeffs_rhs(3);
eqn4 = (226.02*k5 - 122.43*k4 + 226.02*k6) == coeffs_rhs(4);
eqn5 = (0.88*k1 - 4.8*k2 + 0.96*k3 - 171.87) == coeffs_rhs(5);
eqn6 =(0.88 *k4 - 4.8 *k5 + 0.96 *k6)== coeffs_rhs(6);

sol = solve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [k1, k2, k3, k4, k5 ,k6]);
K=double([sol.k1,sol.k2,sol.k3,sol.k4,sol.k5,sol.k6]); % K matrix for pole placement
p=[u1 u2 u3 u4 u5 u6];
    
K_confirm=acker(A,B,p); % confirmation using ackermans
verify_poles=eig(A-B*K);


Q=diag([5000 5000 5000 0 0 0 ]); % Q martix one
Q1=diag([5 50 50 20 700 700]); % Q marix two
R=1; % R matrix

[X,~,~] = care(A,B,Q); % solve the Ricatti Equation

lqr1=inv(R)*B'*X; % LQR controller 1
lqr2=lqr(A,B,Q1,R); % LQR controller 2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Luenberg Full state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms ke1 ke2 ke3 ke4 ke5 ke6
Ke=[ke1; ke2; ke3; ke4; ke5; ke6];

lhs_obs=vpa(det(sI-A+Ke*C));
cl_char_eq_obs=expand((s-u1*5)*(s-u2*5)*(s-u3*5)*(s-u4*5)*(s-u5*5)*(s-u6*5)); % closed loop char eq of Luenberger, 5 times faster

simplified_cl_char_eq_obs=vpa(cl_char_eq_obs);
coeffs_rhs_obs=coeffs(simplified_cl_char_eq_obs);

pretty(collect(lhs_obs));

 eqn1e = 5639.9*ke4 + 503.56*ke5 + 402.85*ke6 == coeffs_rhs_obs(1);
eqn2e = 5639.9*ke1 + 503.56*ke2 + 402.85*ke3 + 5639.9*ke4 + 503.56*ke5 + 402.85*ke6== coeffs_rhs_obs(2);
eqn3e =5639.9 - 5237.0*ke3 - 233.39*ke4 - 182.06*ke5 - 223.8*ke6- 5136.3*ke2 == coeffs_rhs_obs(3);
eqn4e = -233.39 *ke1 - 182.06*ke2 - 223.8*ke3 - 233.39*ke4 - 182.06*ke5 - 223.8*ke6 == coeffs_rhs_obs(4);
eqn5e = 51.331*ke2 + 9.5819*ke3 + ke4 + ke5 + ke6 - 233.39== coeffs_rhs_obs(5);
eqn6e = ke1 + ke2 + ke3 + ke4 + ke5 + ke6 == coeffs_rhs_obs(6);

esol = solve([eqn1e, eqn2e, eqn3e, eqn4e, eqn5e, eqn6e], [ke1, ke2, ke3, ke4, ke5 ,ke6]);

Ke=double([esol.ke1,esol.ke2,esol.ke3,esol.ke4,esol.ke5,esol.ke6]); % observer gain matrix

pe=[u1*5 u2*5 u3*5 u4*5 u5*5 u6*5]; % confirmation observer poles 5 times faster

Ke_confirm=acker(A',C',pe); % confirm observer gain matrix

verify_poles_obs=eig(A-Ke'*C);


%%%%%%%%%%%%%
% For Microphone params
%%%%%%%%%%%%

k=1;
mass=0.5;
damper=200;
R=100;
L=0.0001;
ep=1;
a=0.1;

verify_poles=eig(A-B*K);

%%%%%%%%%%%%%%%

%%%% Microphone params

%%%%%%%%
k=1; % spring const
mass=0.5; %mass
damper=200; % damping coeff
R=100; % resistor val
L=0.0001; % Inductor val
ep=1; % epsilon in air
a=0.1; % area of plates

verify_poles=eig(A-B*K); % verify new poles

%%%%%%%%%%%%%%%%%%%%%%%

% Linearised mic state space

A_=[0 1 0 0; -k/mass -damper/mass (sqrt(ep*2*a*k))/ep*a*mass 0; 0 0 0 1; 0 0 0 -R/L];
B_=[0; 1/mass; 0;0];

% Root locus of double pendulum
%  [num,den] = ss2tf(A,B,C,D);
% 
%  
%  h = tf(num,den)
% rlocus(h) %unstable system

%eigs(A) %poles of the unstable system