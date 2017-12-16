function[feas, xIter,uIter, JIter] = Planner_BatchforMPC_Noa_NoYalmip(u2d,u1d,qx,rou1,rou2,Tx,M,Q, R, N, x0, xL, xU, uL, uU,obstacle,safetyR)
%% Solve Initial guess
% Define state matrix
x = sdpvar(5,N+1);
% Define decision variables
u = sdpvar(1,N);

% Define objective function
objective = x(:,N+1)'*Q*x(:,N+1)+qx*(Tx(M-1+N+1)-x(1,N+1))^2;
for i=1:N
    objective=objective+x(:,i)'*Q*x(:,i)++u(i)'*R*u(i)+qx*(Tx(M-1+i)-x(1,i))^2;
end
    
% Define constraints from 1 to N
constraints = [];
constraints = [constraints x(:,1)==x0];
constraints = [constraints xL<=x(:,N+1)<=xU];
for i = 1:N
    constraints = [constraints  xL<=x(:,i)<=xU uL<=u(i)<=uU ];%upper and lower bound of input and state
    constraints = [constraints  x(:,i+1) == bikeFE_Noa(x(:,i),u(i))];%system dynamics
    if i <= N-1       
%        constraints = [constraints abs(u(2,i+1)-u(2,i))<=u2d abs(u(1,i+1)-u(1,i))<=u1d];
       constraints = [constraints abs(u(i+1)-u(i))<=u2d];
    end
end
% Set options for YALMIP and solver
options = sdpsettings('verbose',0,'solver','IPOPT');
% options = sdpsettings('verbose',0,'solver','fmincon');

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

% judging if the obstacle is seen in a prediction horizon
for xj=xIter{1}(1,:)
    if ((xj<=obstacle(1)+obstacle(3)) && (xj>=obstacle(1)-obstacle(3)))||((xj<=obstacle(4)+obstacle(6)) && (xj>=obstacle(4)-obstacle(6)))
%% Iteration
        k=1;
        tic
        while 1
            x = sdpvar(4,N+1);
            slack =sdpvar(1,N+1);
            % Define decision variables
            u = sdpvar(1,N);

            % Define objective function
%             objective = x(:,N+1)'*Q*x(:,N+1)+qx*(Tx(M-1+N+1)-x(1,N+1))^2+(slack(N+1)^2)*rou1+slack(N+1)*rou2;
%             for i=1:N+1
%                 objective=objective+x(:,i)'*Q*x(:,i)+qx*(Tx(M-1+i)-x(i))^2+(slack(i)^2)*rou1+slack(i)*rou2;
%             end
            objective = x(:,N+1)'*Q*x(:,N+1)+qx*(Tx(M-1+N+1)-x(1,N+1))^2+slack(N+1)*rou2;
            for i=1:N
                objective=objective+x(:,i)'*Q*x(:,i)+qx*(Tx(M-1+i)-x(1,i))^2+u(i)'*R*u(i)+slack(i)*rou2;
            end

            % Define constraints from 1 to N
            constraints = [];
            constraints = [constraints x(:,1)==x0 slack(1)>=0];
            constraints = [constraints xL<=x(:,N+1)<=xU slack(N+1)>=0];
%             l=sqrt(5^2+1.25^2)/2;
%             phi0=atan(1.25/5);
            for i = 1:N
%                 xa=x(1,i)+l*cos(x(4,:)-phi0);
%                 xb=x(1,i)-l*cos(x(4,:)+phi0);
%                 xc=x(1,i)+l*cos(x(4,:)+phi0);
%                 xd=x(1,i)-l*cos(x(4,:)-phi0);
%                 ya=x(2,i)+l*sin(x(4,:)-phi0);
%                 yb=x(2,i)-l*sin(x(4,:)+phi0);
%                 yc=x(2,i)+l*sin(x(4,:)+phi0);
%                 yd=x(2,i)-l*sin(x(4,:)-phi0);
%                 verticies = [xa,ya; xb,yb; xc,yc; xd,yd];
%                 P = Polyhedron('V', verticies);
%                 d=P.distance([obstacle(1);obstacle(2)]);
                constraints = [constraints  xL<=x(:,i)<=xU uL<=u(i)<=uU slack(i)>=0];%upper and lower bound of input and state
                constraints = [constraints  x(:,i+1) == bikeFE_Noa(x(:,i),u(i))];%system dynamics
                constraints = [constraints  (xIter{k}(1,i)-obstacle(1))^2+(xIter{k}(2,i)-obstacle(2))^2+(2*xIter{k}(1,i)-2*obstacle(1))*(x(1,i)-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(2))*(x(2,i)-xIter{k}(2,i))>=(obstacle(3)+safetyR)^2-slack(i)];
                constraints = [constraints  (xIter{k}(1,i)-obstacle(4))^2+(xIter{k}(2,i)-obstacle(5))^2+(2*xIter{k}(1,i)-2*obstacle(4))*(x(1,i)-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(5))*(x(2,i)-xIter{k}(2,i))>=(obstacle(6)+safetyR)^2-slack(i)];
%                 constraints = [constraints  (xIter{k}(1,i)-obstacle(1))*(x(1,i)-obstacle(1))+(xIter{k}(2,i)-obstacle(2))*(x(2,i)-obstacle(2))>=...
%                     0.5*((obstacle(3)+safetyR)^2+(xIter{k}(1,i)-obstacle(1))^2+(xIter{k}(2,i)-obstacle(2))^2)-slack(i)];
%                 constraints = [constraints  (xIter{k}(1,i)-obstacle(1))^2+(xIter{k}(2,i)-obstacle(2))^2+(2*xIter{k}(1,i)-2*obstacle(1))*(d.y(1)-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(2))*(d.y(2)-xIter{k}(2,i))>=(obstacle(3)+safetyR)^2-slack(i)];
%                 constraints = [constraints  (xIter{k}(1,i)-obstacle(1))^2+(xIter{k}(2,i)-obstacle(2))^2+(2*xIter{k}(1,i)-2*obstacle(1))*(xa-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(2))*(ya-xIter{k}(2,i))>=(obstacle(3)+safetyR)^2-slack(i)];
%                 constraints = [constraints  (xIter{k}(1,i)-obstacle(1))^2+(xIter{k}(2,i)-obstacle(2))^2+(2*xIter{k}(1,i)-2*obstacle(1))*(xb-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(2))*(yb-xIter{k}(2,i))>=(obstacle(3)+safetyR)^2-slack(i)];
%                 constraints = [constraints  (xIter{k}(1,i)-obstacle(1))^2+(xIter{k}(2,i)-obstacle(2))^2+(2*xIter{k}(1,i)-2*obstacle(1))*(xc-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(2))*(yc-xIter{k}(2,i))>=(obstacle(3)+safetyR)^2-slack(i)];
%                 constraints = [constraints  (xIter{k}(1,i)-obstacle(1))^2+(xIter{k}(2,i)-obstacle(2))^2+(2*xIter{k}(1,i)-2*obstacle(1))*(xd-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(2))*(yd-xIter{k}(2,i))>=(obstacle(3)+safetyR)^2-slack(i)];
                if i <= N-1       
        %             constraints = [constraints abs(u(2,i+1)-u(2,i))<=0.2 abs(u(1,i+1)-u(1,i))<=0.06];
                   constraints = [constraints abs(u(i+1)-u(i))<=u2d];
                end
            end

            % Set options for YALMIP and solver
            options = sdpsettings('verbose',0,'solver','IPOPT');
%             options = sdpsettings('verbose',0,'solver','fmincon');
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
        toc
        break
        
    else
        continue
    end
end
end