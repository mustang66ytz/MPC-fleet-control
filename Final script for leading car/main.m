%% Main
%Objective is to minimize y, phi, beta and a through start to end where a
%obstacle occurs in the way
clear all 
clear('yalmip')
%QR Tunning
Q=0.1*[0 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
R=[0 0;0 0.1];
rou=1000;
qx=0.1;
%Intial and final position(for batch)
x0=[0;0;5;0];

%Lower and upper constraints of input and state
xL=[0;-4.5;0;-90*pi/180];
xU=[1000;4.5;30;90*pi/180];
uL=-0.6;
uU=0.6;
ud=0.05;
%Obstacle
obstacle=[10;0;1;25;1;1.5;40;0;1;55;-1;1.5];%located at x=25m,y=0m and is a circular obstacle with radius of 2m
safetyR=2;
%% four obstacle_Noa(MPC)
tic
N=5;
M=100;
R=1;
[feas, xOpt, uOpt,JOpt] = Planner_MPC_Noa(u2d,u1d,qx,rou1,rou2,Tx,x0, M, N,Q, R, xL, xU, uL, uU,obstacle,safetyR);
toc
save('four obstacles MPC Test_Noa')