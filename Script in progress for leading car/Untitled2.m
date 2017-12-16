n = 20;
a = 100;
lower = 0;
upper = 2;
temp = exp(linspace(log(1)*a,log(1.05)*a,n));
% re-scale to be between 0 and 1
temp_01 = temp/max(temp) - min(temp)/max(temp);
% re-scale to be between your limits (i.e. 1 and 1.05)
out = temp_01*(upper-lower) + lower;
plot(diff(out),diff(out),'o')
%%
a=logspace(0,log10(101),10)-1;
a=a/100;
b=logspace(-log10(101),0,10)-1;
b=b/100;
x=linspace(0,0,10);
plot(b,x,'-o')
%%
clear
a=linspace(0,pi/2,10);
b=sin(a);
c=linspace(pi/2,0,10);
d=cos(a)+1;
% plot(a,b,'-o')
x=linspace(0,0,20);
haha=cat(2,b,d);

plot(haha,x,'-o')



