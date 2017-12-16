function[feas, xIter,uIter, JIter] = Reformulation_1_BatchMPC(P0,R, N, x0, xL, xU, uL, uU,obstacle,safetyR)
%% Solve Initial guess
% Define state matrix
x = sdpvar(5,N+1);% x,y,v,phi,slack
xbar=sdpvar(5,5,N+1);%semidefinite relaxation variable

% Define decision variables
u = sdpvar(2,N);
ubar=sdpvar(2,2,N);%semidefinite relaxation variable
% Tracking 
Tx=linspace(0,50,N+1);

% Define objective function
rou=1000;
q0=[-2*Tx(N+1);0;0;0;rou];
% objective = x(:,N+1)'*P0*x(:,N+1)+q0'*x(:,N+1)+Tx(N+1)^2;
objective = trace(xbar(:,:,N+1)*P0)+q0'*x(:,N+1)+Tx(N+1)^2;
for i=1:N
    q0=[-2*Tx(i);0;0;0;rou];
    objective=objective+trace(xbar(:,:,i)*P0)+trace(ubar(:,:,i)*R)+q0'*x(:,i)+Tx(i)^2;
end
    
% Define constraints from 1 to N
constraints = [];
constraints = [constraints x(1:4,1)==x0 x(5,1)>=0];
constraints = [constraints xL<=x(1:4,N+1)<=xU x(5,N+1)>=0];
xcat=cat(2,cat(1,xbar(:,:,N+1),x(:,N+1)'),cat(1,x(:,N+1),[1]));
constraints = [constraints  eig(xcat)>=0];%SDP loose constraints
% constraints = [constraints  xcat>=eye(6) ];%SDP loose constraints

for i = 1:N
    constraints = [constraints  xL<=x(1:4,i)<=xU uL<=u(:,i)<=uU x(5,i)>=0];%upper and lower bound of input and state
    constraints = [constraints  x(1:4,i+1) == bikeFE(x(1:4,i),u(:,i))];%system dynamics
    xcat=cat(2,cat(1,xbar(:,:,i),x(:,i)'),cat(1,x(:,i),[1]));
    ucat=cat(2,cat(1,ubar(:,:,i),u(:,i)'),cat(1,u(:,i),[1]));
    constraints = [constraints  eig(xcat)>=0 ];
    constraints = [constraints  eig(ucat)>=0 ];
%     xe=eig(xcat);
%     ue=eig(ucat);
%     constraints = [constraints  xe(1)>0 xe(2)>0 xe(3)>0 xe(4)>0 xe(5)>0 xe(6)>0];%SDP loose constraints
%     constraints = [constraints  ue(1)>0 ue(2)>0];%SDP loose constraints
   
    if i <= N-1       
        constraints = [constraints abs(u(2,i+1)-u(2,i))<=0.2 abs(u(1,i+1)-u(1,i))<=0.06];
    end
end
% Set options for YALMIP and solver
%options = sdpsettings('verbose',0,'solver','IPOPT');
options = sdpsettings('verbose',0,'solver','quadprog');
optimize(constraints, objective, options);

%diagnostics
diagnostics=optimize(constraints, objective, options);
if diagnostics.problem == 0 
    xIter{1}=value(x);
    uIter{1}=value(u); 
    JIter{1}=value(objective);
    feas=true;
    fprintf('Intial is solved with Jopt being %f \n',JIter{1})
else
    feas=false;
    xIter{1}=[];
    uIter{1}=[];
    JIter{1}=[];
    return
end

%% Iteration
k=1;
while 1
    % Define state matrix
    x = sdpvar(5,N+1);% x,y,v,phi,slack
    xbar=sdpvar(5,5,N+1);%semidefinite relaxation variable
    % Define decision variables
    u = sdpvar(2,N);
    ubar=sdpvar(2,2,N);%semidefinite relaxation variable
    % Tracking 
    Tx=linspace(0,50,N+1);

    % Define objective function
    rou=1000;
    P0=[1 0 0 0 0;0 500 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0];
    q0=[-2*Tx(N+1);0;0;0;rou];
    objective = trace(xbar(:,:,N+1)*P0)+q0'*x(:,N+1)+Tx(N+1)^2;
%     objective = x(:,N+1)'*P0*x(:,N+1)+q0'*x(:,N+1)+Tx(N+1)^2;
    for i=1:N
        q0=[-2*Tx(i);0;0;0;rou];
        objective=objective+trace(xbar(:,:,i)*P0)+trace(ubar(:,:,i)*R)+q0'*x(:,i)+Tx(i)^2;
    end

    % Define constraints from 1 to N
    constraints = [];
    constraints = [constraints x(1:4,1)==x0 x(5,1)>=0];
    constraints = [constraints xL<=x(1:4,N+1)<=xU x(5,N+1)>=0];
    xcat=cat(2,cat(1,xbar(:,:,N+1),x(:,N+1)'),cat(1,x(:,N+1),[1]));
%     constraints = [constraints  xcat>=eye(6) ];%SDP loose constraints
    constraints = [constraints  eig(xcat)>=0];
    for i = 1:N
        constraints = [constraints  xL<=x(1:4,i)<=xU uL<=u(:,i)<=uU x(5,i)>=0];%upper and lower bound of input and state
        constraints = [constraints  x(1:4,i+1) == bikeFE(x(1:4,i),u(:,i))];%system dynamics
        constraints = [constraints  (xIter{k}(1,i)-obstacle(1))^2+(xIter{k}(2,i)-obstacle(2))^2+(2*xIter{k}(1,i)-2*obstacle(1))*(x(1,i)-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(2))*(x(2,i)-xIter{k}(2,i))>=(obstacle(3)+safetyR)^2-slack(i)];
        xcat=cat(2,cat(1,xbar(:,:,i),x(:,i)'),cat(1,x(:,i),[1]));
        ucat=cat(2,cat(1,ubar(:,:,i),u(:,i)'),cat(1,u(:,i),[1]));
        constraints = [constraints  eig(xcat)>=0];
        constraints = [constraints  eig(ucat)>=0];
        if i <= N-1       
            constraints = [constraints abs(u(2,i+1)-u(2,i))<=0.2 abs(u(1,i+1)-u(1,i))<=0.06];
        end
    end
    % Set options for YALMIP and solver
    options = sdpsettings('verbose',0,'solver','quadprog');
  %   options = sdpsettings('verbose',0,'solver','sedumi');

    % Solve
    optimize(constraints, objective, options);

    %diagnostics
    diagnostics=optimize(constraints, objective, options);
    if diagnostics.problem == 0 
        xIter{k+1}=value(x);
        uIter{k+1}=value(u); 
        JIter{k+1}=value(objective);
        feas=true;
        fprintf('%d th iteration is solved with Jopt being %f \n',k,JIter{k+1})
    else
        feas=false;
        fprintf('%d th iteration is infeasible \n',k)
%         xIter{k+1}=[];
%         uIter{k+1}=[];
%         JIter{k+1}=[];
        return
    end
    %End of iteration
    if abs(JIter{k+1}-JIter{k})<1
        disp('Iterations converge at epsilon=1')
        break
    end
    if k>25
        disp('Too many iteration(More than25)! Abandon iterations and jump to next timestep \n')
        break
    end
    k=k+1;
end
end