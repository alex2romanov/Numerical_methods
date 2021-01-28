N = 10;
mu = 3;
X = zeros(1, N+1);
Y = zeros(1, N+1);
L = zeros(1, N+1);
% после 5000 итерации получил формулу для аргумента 
for w = 1:N+1
        z = ((2+3*mu)*exp(-i*w)+8+(2-3*mu)*exp(i*w))/((2-3*mu)*exp(-i*w)+8+(2+3*mu)*exp(i*w))
        X(w) = angle(z)
        
        Y(w) = -2*atan(3*mu*sin(w)/(2*cos(w)+mu*mu*cos(w)-mu*mu+4))
        %L(w) = atan(-4*mu*sin(w)/(4-mu*mu*sin(w)*sin(w)))-pi
end
