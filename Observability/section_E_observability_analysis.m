%%
% Section E Code: Study of Observability: 
syms M m1 m2 l1 l2 g;
%Declaring variables;
M= 1000; %Crane mass
m1= 100; % Load 1 mass
m2= 100; % Load 2 mass
l1= 20; % Cable length of Load 1
l2= 10; % Cable length of Load 2
g= 9.81;
% Creating a linearised state space equation using A and B matrices
A=[0 1     0               0       0            0;
   0 0 -(m1*g)/M           0 -(m2*g)/M          0;
   0 0     0               1       0            0;
   0 0 -((M+m1)*g)/(M*l1)  0 -(m2*g)/(M*l1)     0;
   0 0     0               0       0            1;
   0 0 -(m1*g)/(M*l2)      0 -(g*(M+m2))/(M*l2) 0];
B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)]; %Initializing the Linearized B matrix
C1 = [1 0 0 0 0 0]; %Case 1: Only x is observed
C2 = [0 0 1 0 0 0; 0 0 0 0 1 0]; %Case 2: theta1 and theta2 are observed.
C3 = [1 0 0 0 0 0; 0 0 0 0 1 0]; % Case 3: x and theta2 are observed. 
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]; %Case 4: x,theta1 and theta2 are observed. 
%Individual matrix to check each Observability Condition
Observability1 = obsv(A,C1);
Observability2 = obsv(A,C2);
Observability3 = obsv(A,C3);
Observability4 = obsv(A,C4);
rankArray = [rank(Observability1),rank(Observability2),rank(Observability3),rank(Observability4)];
disp("Matrix A of Linearized model is: "); disp(A);
disp("Matrix B of Linearized model is: "); disp(B);
disp("Rank of each Observability matrix");
disp(rankArray);
