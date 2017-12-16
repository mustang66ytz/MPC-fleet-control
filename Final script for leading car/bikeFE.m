function zp = bikeFE(z, u)
    % z(1) = x
    % z(2) = y
    % z(3) = v
    % z(4) = phi
    % u = steering
    % acceleration is constant zero 
    TS = 0.2;
    lf = 3;
    lr = 3;
    beta = atan(lr*tan(u)/(lf+lr));
    zp = sdpvar(4,1);
    
    zp(1) = z(1)+TS*z(3)*cos(z(4)+beta);
    zp(2) = z(2)+TS*z(3)*sin(z(4)+beta);
    zp(3) = z(3);
    zp(4) = z(4)+TS*z(3)*sin(beta)/lr;
    
end