%%
% Section D : Checking that the system is controllable after value substitution

% Declaring variables;
M= 1000; %Crane mass
m1= 100; % Load 1 mass
m2= 100; % Load 2 mass
l1= 20; % Cable length of Load 1
l2= 10; % Cable length of Load 2
g= 9.81;

A=[0 1     0               0       0            0;
   0 0 -(m1*g)/M           0 -(m2*g)/M          0;
   0 0     0               1       0            0;
   0 0 -((M+m1)*g)/(M*l1)  0 -(m2*g)/(M*l1)     0;
   0 0     0               0       0            1;
   0 0 -(m1*g)/(M*l2)      0 -(g*(M+m2))/(M*l2) 0];
disp("A matrix after values substitution: "); disp(A);

% Matrix B
B=[0; 1/M; 0; 1/(M*(l1)); 0; 1/(M*l2)];
disp("B matrix after value substitution : ");disp(B)

disp("Rank of controllability matrix: "); 
disp(rank(ctrb(A,B)));

disp("Controllability Matrix is = ")
disp(ctrb(A,B));
disp(det(ctrb(A,B)));

disp("Eigen Values of A");
disp(eig(A));

%%
% Section to check the initial state of the system without LQR. 
% The initial conditions are as follows.
X_0 = [0;0;10;0;20;0];
% We assume the values of Q and R.
Q=[ 1 0   0  0  0  0;
    0 10  0  0  0  0;
    0 0 1000 0  0  0;
    0 0   0 10 0   0;
    0 0   0  0 1000 0;
    0 0   0  0  0 10];
R=0.001;
% Q and R are a part of the cost function of LQR controller
% It is an trade-off between Q and R so we use them both to
% develop a system as per our priorities. 
% Assumption: C matrix is a direct representation of the output
% matrix, which makes D=0
disp("Using MATLAB inbuild function ss to check the initial response of the system.");
disp("Using Q and R matrix of our assumption:");
C = eye(6); D = 0;
sys1 = ss(A,B,C,D);
% ss is the MATLAB function for
% calculating the state space representation of the system
disp("Starting grid"); 
figure
initial(sys1,X_0);

grid on


%%
% Design of the LQR controller is done in this section. 
% Q and R are a part of the cost function of LQR controller
% It is an trade-off between Q and R so we use them both to
% develop a system as per our priorities. 
% Assumption: C matrix is a direct representation of the output
% matrix, which makes D=0
disp("Using MATLAB inbuild function ss to check the initial response of the system.");
disp("Using Q and R matrix of our assumption:");
C = eye(6); D = 0; 
sys1 = ss(A,B,C,D); 
% ss is the MATLAB function for
% calculating the state space representation of the system
disp("Starting grid"); 
figure
initial(sys1,X_0);

grid on

%In-built MATLAB code for LQR Controllers
[K_Gain_mat, Po_def_mat, Poles] = lqr(A,B,Q,R);
disp("K gain matrix using lqr function of MATLAB: ");
disp(K_Gain_mat);

disp("Solution of the associated Ricatti equation using lqr function of MATLAB: ");
disp(Po_def_mat);

disp("Pole placement using lqr function of MATLAB: ");
disp(Poles);

disp("Using MATLAB inbuild function ss to check the initial response of the system.");
disp("Using Q and R matrix calculated by lqr matrix.");

sys2 = ss(A-(B*K_Gain_mat),B,C,D);
%Using the K matrix to define ss
figure
initial(sys2,X_0)
grid on
