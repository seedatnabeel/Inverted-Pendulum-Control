clc
clear all

%%%%%%%%%%%%%%%

%%%%%   INIT PENDULUM PARAMETERS

M = 1; %mass of cart
m = 0.2; %mass of pendulum
l = 0.3; %length of pendulum
b = 0.1; % damping
g=9.81; % gravity
F=1; % force

%%%%%%%%

% LINEARISED PENDULUM MATRIX
A= [0 1 0 0; (g*(M+m))/(l*M) 0 0 b/(l*M); 0 0 0 1; -m*g/M 0 0 -b/M];

B = [0;-1/(l*M); 0 ; 1/M];

C=[1 1 1 1];
D=0;

unstable_poles=eigs(A); %calculate unstable poles
%[num,den] = ss2tf(A,B,C,D);
%
%   h = tf(num,den);
%   figure
%   title('Linear Single Link Pendulum Root Locus');
%  rlocus(h)
%

%%% prove controllability
single_ctrb_matrix=horzcat(B, A*B, A^2*B, A^3*B);
unctrb_states=length(A) - rank(single_ctrb_matrix); %test controllability

%confirm with ctrb
ctrb_func= ctrb(A,B);

%%%prove observavbility
single_obs_matrix=obsv(A,C);


%%%%%%%%%%%%%%%

% Pole placement

%%%%%%%%%%%%%%%
syms s k1 k2 k3 k4 k5 k6;
sI=[s 0 0 0 ; 0 s 0 0  ; 0 0 s 0 ; 0 0 0 s];

desired_roots_non_repeated=roots([1,2*0.59115*3.383,3.383*3.383]); %find non repeated roots
desired_roots_repeated=roots([1,2*0.59115*3.383,3.383*3.383]); %find repeated roots

% Pole placement matrix algebra
u1=desired_roots_non_repeated(1);
u2=desired_roots_non_repeated(2);
u3=desired_roots_non_repeated(1)*5;
u4=desired_roots_non_repeated(2)*5;
ctlr_poles=[u1 u2 u3 u4];
digits(5);
cl_char_eq=expand((s-u1)*(s-u2)*(s-u3)*(s-u4)); % desired closed loop eq
simplified_cl_char_eq=vpa(cl_char_eq);
K=[k1 k2 k3 k4];
lhs=vpa(det(sI-A+B*K));% calc determinant

pretty(simplified_cl_char_eq);
pretty(collect(lhs));

coeffs_rhs=coeffs(simplified_cl_char_eq);

eqn1 = - 32.7 *k3 == coeffs_rhs(1);
eqn2 = (- 32.7 *k4 - 3.27)== coeffs_rhs(2);
eqn3 = (k3 - 3.3333 *k1 - 39.24) == coeffs_rhs(3);
eqn4 = (k4 - 3.3333 *k2 + 0.1)  == coeffs_rhs(4);

sol = solve([eqn1, eqn2, eqn3, eqn4], [k1, k2, k3, k4]);

K_non_rep=double([sol.k1,sol.k2,sol.k3,sol.k4]); % Non repeated poles k matrix
KA=K_non_rep;
p2=[desired_roots_repeated(1),desired_roots_repeated(2) desired_roots_repeated(1) desired_roots_repeated(2) ];

K_rep=acker(A,B,p2); % repeated poles k matrix
KB=K_non_rep;

%%%%% Different tested Q and R parameters
Q1=diag([500 0 5000 0 ]);
Q2=diag([5000 0 100 0 ]);

R1=1;
R2=2;

%%%%%% 4 Different
K1=lqr(A,B,Q1,R1);
K2=lqr(A,B,Q2,R1);
K3=lqr(A,B,Q1,R2);
K4=lqr(A,B,Q2,R2);

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

