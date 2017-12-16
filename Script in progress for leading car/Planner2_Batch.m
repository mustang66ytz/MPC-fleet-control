function[feas, xIter,uIter, JIter] = Planner2_Batch(Q, R, N, x0,xN, xL, xU, uL, uU,obstacle,safetyR)
%% Solve Initial guess
% Define state matrix
x = sdpvar(4,N+1);
slack =sdpvar(1,N+1);
rou=1000;
% Define decision variables
u = sdpvar(2,N);

% Tracking 
Tx=linspace(0,50,N+1);

% Define objective function
objective = x(:,N+1)'*Q*x(:,N+1)+(Tx(N+1)-x(1,N+1))^2+slack(N+1)*rou;
for i=1:N
    objective=objective+x(:,i)'*Q*x(:,i)++u(:,i)'*R*u(:,i)+(Tx(i)-x(i))^2+slack(i)*rou;
end
    
% Define constraints from 1 to N
constraints = [];
constraints = [constraints x(:,1)==x0 slack(1)>=0];
constraints = [constraints x(:,N+1)==xN slack(N+1)>=0];
for i = 1:N
    constraints = [constraints  xL<=x(:,i)<=xU uL<=u(:,i)<=uU slack(i)>=0];%upper and lower bound of input and state
    constraints = [constraints  x(:,i+1) == bikeFE(x(:,i),u(:,i))];%system dynamics
%     constraints = [constraints  norm([x(1,i+1)-x(1,i) x(2,i+1)-x(2,i)])<=norm([x(1,i)-obstacle(1) x(2,i)-obstacle(2)])-safetyR-obstacle(3)];%
    if i <= N-1       
        constraints = [constraints abs(u(2,i+1)-u(2,i))<=0.2 abs(u(1,i+1)-u(1,i))<=0.06];
    end
end
% Set options for YALMIP and solver
options = sdpsettings('verbose',0,'solver','IPOPT');

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
    x = sdpvar(4,N+1);
    slack =sdpvar(1,N+1);
    rou=1000;
    % Define decision variables
    u = sdpvar(2,N);

    % Tracking 
    Tx=linspace(0,50,N+1);

    % Define objective function
    objective = x(:,N+1)'*Q*x(:,N+1)++u(:,i)'*R*u(:,i)+(Tx(N+1)-x(1,N+1))^2+slack(N+1)*rou;
    for i=1:N
        objective=objective+x(:,i)'*Q*x(:,i)+(Tx(i)-x(i))^2+slack(i)*rou;
    end

    % Define constraints from 1 to N
    constraints = [];
    constraints = [constraints x(:,1)==x0 slack(1)>=0];
    constraints = [constraints x(:,N+1)==xN slack(N+1)>=0];
    for i = 1:N
        constraints = [constraints  xL<=x(:,i)<=xU uL<=u(:,i)<=uU slack(i)>=0];%upper and lower bound of input and state
        constraints = [constraints  x(:,i+1) == bikeFE(x(:,i),u(:,i))];%system dynamics
        constraints = [constraints  (xIter{k}(1,i)-obstacle(1))^2+(xIter{k}(2,i)-obstacle(2))^2+(2*xIter{k}(1,i)-2*obstacle(1))*(x(1,i)-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(2))*(x(2,i)-xIter{k}(2,i))>=(obstacle(3)+safetyR)^2-slack(i)];
        constraints = [constraints  (xIter{k}(1,i)-obstacle(4))^2+(xIter{k}(2,i)-obstacle(5))^2+(2*xIter{k}(1,i)-2*obstacle(4))*(x(1,i)-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(5))*(x(2,i)-xIter{k}(2,i))>=(obstacle(6)+safetyR)^2-slack(i)];
        if i <= N-1       
            constraints = [constraints abs(u(2,i+1)-u(2,i))<=0.2 abs(u(1,i+1)-u(1,i))<=0.06];
        end
    end

    % Set options for YALMIP and solver
    options = sdpsettings('verbose',0,'solver','IPOPT');

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
    k=k+1;
end
end