syms a0 b0 a1 b1 p0 p1 q0 q1 V h t;

p0 = -1;

A = 2*a0+b0+2*a1+b1;

B = (2*a0+b0)*t-2*p0-q0-2*p1-q1;

C = (2*a0+2*a1)*h^2+12*V*(2*p0+2*p1);

D = (2*a0+2*a1)*h^2+2*V*(2*p0+q0+2*p1+q1);

E = 2*a0*h^2*t-(2*p0+2*p1)*h^2+2*V*t*(2*p0+q0);

F = (2*a0+b0)*t-2*(2*p0+q0);

G = 2*a0*h^2*t^2-4*p0*h^2*t+2*V*(2*p0+q0)*t^2;


s  = solve(A,B, C, D, E, F, G, a1, b1,  a0 ,b0,p1,q0, q1);


a1 = s.a1;
b1 = s.b1;
a0 = s.a0;
b0 = s.b0;
p1 = s.p1;
q0 = s.q0;
q1 = s.q1;
X = linspace(0, 2*pi, n1+1);
X = X(1:n1)';
u0 = sin(X)*cos(0);
plot(X, u0);
