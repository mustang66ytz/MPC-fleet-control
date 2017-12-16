figure
plot(xOpt(1,1:60),uOpt(:),'-o')
ylim([-0.8 0.8])
hline1=refline(0,0.6);
hline2=refline(0,-0.6);
hline1.LineStyle=('--');
hline2.LineStyle=('--');
xlabel('x(m)')
ylabel('steering angle(radian)')
%%
figure
N=[5,8,11];
T=[1.9,4.2,9.3];
plot(N,T,'-o')