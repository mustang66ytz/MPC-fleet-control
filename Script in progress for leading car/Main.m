%% Main
%Objective is to minimize y, phi, beta and a through start to end where a
%obstacle occurs in the way
clear all 
clear('yalmip')
%QR Tunning
Q=0.1*[0 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
R=[0 0;0 0.1];
rou1=100;
rou2=1000;
qx=0.1;
%Intial and final position(for batch)
x0=[0;0;5;0];
xN=[100;0;5;0];%(for batch)

%Lower and upper constraints of input and state
xL=[0;-4.5;0;-90*pi/180];
xU=[1000;4.5;30;90*pi/180];
uL=[-5.4;-0.6];
uU=[2.7;0.6];
u2d=0.05;
u1d=0.5;
%Obstacle
obstacle=[10;0;1;25;1;1.5;40;0;1;55;-1;1.5];%located at x=25m,y=0m and is a circular obstacle with radius of 2m
safetyR=2;

% Tracking 
% T1=linspace(0,22,23);
% T2=linspace(23,27,10);
% T3=linspace(28,72,45);
% T4=linspace(73,77,10);
% T5=linspace(78,1000,923);
% Tx=cat(2,T1,T2,T3,T4,T5);
Tx=linspace(0,1000,1001);
%% One obstacle(Batch)
% Computation
tic
N=50;
[feas, xIter,uIter,JIter]  = Planner_Batch(Q, R, N, x0,xN, xL, xU, uL, uU,obstacle,safetyR);
toc
% Plot
figure
plot(xIter{end}(1,:),xIter{end}(2,:),'-o')
axis equal
viscircles([obstacle(1) obstacle(2)],safetyR+obstacle(3))
viscircles([obstacle(1) obstacle(2)],obstacle(3))
hold on
plot(xIter{end}(1,:),xIter{end}(3,:))
plot(xIter{end}(1,:),xIter{end}(4,:))
plot(xIter{end}(1,1:N),uIter{end}(1,:))
plot(xIter{end}(1,1:N),uIter{end}(2,:))
legend('y','v','phi','a','beta')
% Save
save('One obstacle Test')
clear all
%% Two obstacles(Batch)
% Computation
tic
N=50;
[feas, xIter,uIter,JIter]  = Planner2_Batch(Q, R, N, x0,xN, xL, xU, uL, uU,obstacle,safetyR);
toc
% Plot
figure
plot(xIter{end}(1,:),xIter{end}(2,:),'-o')
axis equal
viscircles([obstacle(1) obstacle(2)],safetyR+obstacle(3))
viscircles([obstacle(1) obstacle(2)],obstacle(3))
viscircles([obstacle(4) obstacle(5)],safetyR+obstacle(6))
viscircles([obstacle(4) obstacle(5)],obstacle(6))
hold on
plot(xIter{end}(1,:),xIter{end}(3,:))
plot(xIter{end}(1,:),xIter{end}(4,:))
plot(xIter{end}(1,1:N),uIter{end}(1,:))
plot(xIter{end}(1,1:N),uIter{end}(2,:))
legend('y','v','phi','a','beta')
% Save
save('Two obstacles Test')
clear all
%% One obstacle(MPC)
tic
N=5;
M=30;
[feas, xOpt, uOpt,JOpt] = Planner_MPC(u2d,u1d,qx,rou1,rou2,Tx,x0, M, N,Q, R, xL, xU, uL, uU,obstacle,safetyR);
toc
save('One obstacle MPC Test')
%% One obstacle_Noa(MPC)
tic
N=5;
M=50;
R=1;
[feas, xOpt, uOpt,JOpt] = Planner_MPC_Noa(u2d,u1d,qx,rou1,rou2,Tx,x0, M, N,Q, R, xL, xU, uL, uU,obstacle,safetyR);
toc
save('One obstacle MPC Test_Noa')
%% four obstacle_Noa(MPC)
tic
N=5;
M=100;
R=1;
[feas, xOpt, uOpt,JOpt] = Planner_MPC_Noa(u2d,u1d,qx,rou1,rou2,Tx,x0, M, N,Q, R, xL, xU, uL, uU,obstacle,safetyR);
toc
save('four obstacles MPC Test_Noa')
%% Two obstacles(MPC)
tic
N=10;
M=50;
% xL=[0;-5;0;-90*pi/180];
% xU=[1000;5;10;90*pi/180];
% uL=[-5.4;-0.6];
% uU=[2.7;0.6];
[feas, xOpt, uOpt,JOpt] = Planner2_MPC(x0, M, N,Q, R, xL, xU, uL, uU,obstacle,safetyR);
toc
save('Two obstacles MPC Test')
clear all
%% QCQPS formulation of One obstacle(MPC)
P0=[1 0 0 0 0;0 500 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0];
tic
N=10;
M=50;
xU=[1000;5;10;90*pi/180];
[feas, xOpt, uOpt,JOpt] = Reformulation_1_MPC(x0, M, N,P0, R, xL, xU, uL, uU,obstacle,safetyR);
toc
save('One obstacle MPC by SDPR')

%% One obstacle considering shape(Batch)
% tic
% N=50;
% [feas, xIter,uIter, JIter] = Planner_Batch_shape(Q, R, N, x0,xN, xL, xU, uL, uU,obstacle,safetyR);
% toc
% % Plot
% figure
% plot(xIter{end}(1,:),xIter{end}(2,:),'-o')
% axis equal
% viscircles([obstacle(1) obstacle(2)],safetyR+obstacle(3))
% viscircles([obstacle(1) obstacle(2)],obstacle(3))
% hold on
% plot(xIter{end}(1,:),xIter{end}(3,:))
% plot(xIter{end}(1,:),xIter{end}(4,:))
% plot(xIter{end}(1,1:N),uIter{end}(1,:))
% plot(xIter{end}(1,1:N),uIter{end}(2,:))
% legend('y','v','phi','a','beta')
% % Save
% save('One obstacle Test')
% clear all