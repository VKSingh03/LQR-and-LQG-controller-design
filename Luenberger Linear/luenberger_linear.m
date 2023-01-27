M=1000;%Mass of the Crane
m_1=100;%mass of Load 1
m_2=100;%mass of Load 2
l_1=20;%length of the string of Load 1
l_2=10;%length of the string of Load 2
g=9.81; %declaring the value of the accelertaion due to gravity in m/s^2
%Defining our matrices as follows
A=[0 1 0 0 0 0; 
    0 0 -(m_1*g)/M 0 -(m_2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m_1)*g)/(M*l_1) 0 -(m_2*g)/(M*l_1) 0;
    0 0 0 0 0 1;
    0 0 -(m_1*g)/(M*l_2) 0 -(g*(M+m_2))/(M*l_2) 0];
B=[0; 1/M; 0; 1/(M*l_1); 0; 1/(M*l_2)];
% From previous case, we have determined that only C1, C3 and C4 were
% observable. Hence, we are going to consider only those 3 cases.
C_1 = [1 0 0 0 0 0];  %Corresponding to x component
C_3 = [1 0 0 0 0 0; 0 0 0 0 1 0]; %cooresponding to x and theat2
C_4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]; %cooresponding to x, theta1 and theat2
D = 0;%declarring the D matrix to be zero
% Cosidering the same Q and R matrices chosen before in our code
Q=[1 0 0 0 0 0;
   0 10 0 0 0 0;
   0 0 1000 0 0 0;
   0 0 0 10 0 0;
   0 0 0 0 1000 0;
   0 0 0 0 0 10];
R=0.01; %%these are the cost variables from LQR
% Initial Conditions for Leunberger observer - 12 state variables, 
% 6 actual + 6 estimates
x0=[0,0,30,0,60,0,0,0,0,0,0,0];
% state variables order = [x,dx,theta_1,dthetsa_1,theta_2,dtheta_2,
% estimates taken in the same order]
% For pole placement, lets choose eigen values with negative real part
poles=[-1;-2;-3;-4;-5;-6];
% Calling LQR function to obtain K matrix
K=lqr(A,B,Q,R);
% Framing L for all three cases where output is observable
% Using the pole placement funciton built into MATLAB
L_1 = place(A',C_1',poles)' %L_1 should be a 6x1 matrix
L_3 = place(A',C_3',poles)' %L_3 should be a 6x2 matrix
L_4 = place(A',C_4',poles)' %L_4 should be a 6x3 matrix

A_c1 = [(A-B*K) B*K; % Luenberger A matrix
        zeros(size(A)) (A-L_1*C_1)];
B_c = [B;zeros(size(B))];% Luenberger B matrix
C_c1 = [C_1 zeros(size(C_1))];% Luenberger C matrix

A_c3 = [(A-B*K) B*K;% Luenberger A matrix
        zeros(size(A)) (A-L_3*C_3)];
C_c3 = [C_3 zeros(size(C_3))];% Luenberger C matrix

A_c4 = [(A-B*K) B*K;% Luenberger A matrix
        zeros(size(A)) (A-L_4*C_4)];
C_c4 = [C_4 zeros(size(C_4))];% Luenberger C matrix

sys_1 = ss(A_c1, B_c, C_c1,D);%%MATLAB function to output statespace equations
figure % to launch a new figure WINDOW
initial(sys_1,x0)%MATLAB function to check the initial response of the system
figure
step(sys_1)%Gives the step response output

sys_3 = ss(A_c3, B_c, C_c3,D);%MATLAB function to output statespace equations
figure 
initial(sys_3,x0)%MATLAB inbuilt function to check the initial response of the system
figure 
step(sys_3)%Gives the step response output

sys_4 = ss(A_c4, B_c, C_c4, D);%%MATLAB function to output statespace equations
figure 
initial(sys_4,x0)%MATLAB inbuilt function to check the initial response of the system
figure 
step(sys_4)%Gives the step response output

%grid on
