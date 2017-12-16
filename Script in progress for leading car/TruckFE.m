%% ME C231A HW 3: Optimization Part 2 and Finite-Time Optimal Control
% Edward Zhu
% 9/29/2016

function zp = TruckFE(z, u)
    % z(1) = x
    % z(2) = y
    % z(3) = v
    % z(4) = phi
    % u(1) = a
    % u(2) = beta

    TS = 0.2;
    lf = 5;
    lr = 5;
    %beta = atan2(lr*tan(u(2)),(lf+lr));
    %zp = zeros(4,1);
    zp = sdpvar(4,1);
    
    zp(1) = z(1)+TS*z(3)*cos(z(4)+u(2));
    zp(2) = z(2)+TS*z(3)*sin(z(4)+u(2));
    zp(3) = z(3)+TS*u(1);
    zp(4) = z(4)+TS*z(3)*sin(u(2))/lr;
    
%     xp = x + TS*v*cos(psi+beta);
%     yp = y + TS*v*sin(psi+beta);
%     vp = v + TS*a;
%     psip = psi + TS*v*sin(beta)/lr;
end