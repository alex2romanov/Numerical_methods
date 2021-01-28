syms mu;
mu = 1;
w = linspace(0, 5, 500);
l = mu*(-w);
hold on
plot(w, l)

z = (mu*exp(-j*w)-mu*exp(j*w)+4)/(mu*exp(j*w)-mu*exp(-j*w)+4);
N= 500;
X = zeros(1, N+1);
for k = 1:N+1
    X(k)=(mu*exp(-j*(k-1)/100)-mu*exp(j*(k-1)/100)+4)/(mu*exp(j*(k-1)/100)-mu*exp(-j*(k-1)/100)+4);
end
Y = real(X)

Y = Y(1:500)

plot(w, Y)
xlabel ('w')
ylabel ('arg')
