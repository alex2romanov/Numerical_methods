syms w mu
mu = 100
f = -2*atan(sin(mu*w/2))+mu*w;
%T6 = taylor(f, w);
Y_kn = zeros(1, 4);
Y_comp = zeros(1, 4);
w1 = [0.0001, 0.001, 0.01, 0.1]
for i = 1: 4
    Y_kn(i) =  -2*atan(mu*sin(w1(i))/2)+mu*w1(i)
    Y_comp(i) = -2*atan(3*mu*sin(w1(i))/(2*cos(w1(i))+mu*mu*cos(w1(i))-mu*mu+4))+mu*w1(i)
end

%T8 = taylor(f, x, 'Order', 8);

%g= -2*atan(3*t*sin(x)/(2*cos(x)+t^2*cos(x)-t^2+4))+t*x
%T7 = taylor(g, x);

%plot(log(T6), log(w))
hold on
loglog(w1, Y_kn)
loglog(w1, Y_comp)
hold off