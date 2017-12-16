%% ME C231A HW 3: Optimization Part 2 and Finite-Time Optimal Control
% Edward Zhu
% 9/29/2016

function zp = bikeFE_Noa(z, u)
    % z(1) = x
    % z(2) = y
    % z(3) = v
    % z(4) = phi
    % u(1) = a = constant
    % u(2) = steering

    TS = 0.2;
    lf = 3;
    lr = 3;
    u=1.5*u;
    beta = atan(lr*tan(u)/(lf+lr));
    %zp = zeros(4,1);
    zp = sdpvar(4,1);
    
    
    zp(1) = z(1)+TS*z(3)*cos(z(4)+beta);
    zp(2) = z(2)+TS*z(3)*sin(z(4)+beta);
%     zp(3) = z(3)+TS*u(1);
    zp(3) = z(3);
    zp(4) = z(4)+TS*z(3)*sin(beta)/lr;
    
%     xp = x + TS*v*cos(psi+beta);
%     yp = y + TS*v*sin(psi+beta);
%     vp = v + TS*a;
%     psip = psi + TS*v*sin(beta)/lr;
end