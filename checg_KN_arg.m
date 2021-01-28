N = 10;
mu = 3;
X = zeros(1, N+1);
Y = zeros(1, N+1);
L = zeros(1, N+1);
for w = 1:N+1
        z = (mu*exp(-i*w)-mu*exp(i*w)+4)/(mu*exp(i*w)-mu*exp(-i*w)+4)
        X(w) = angle(z)
        
        Y(w) = -2*atan(mu*sin(w)/2)
        L(w) = atan(-4*mu*sin(w)/(4-mu*mu*sin(w)*sin(w)))-pi
end

% check - success