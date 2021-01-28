
syms a1 b1 c1 a0 b0 c0 p0 p1 q0 q1 r0 r1  h t V w z;

p0 = 1;

A = a1+b1+c1+a0+b0+c0;

B = (a1+b1+c1)*t-p1-q1-r1-p0-q0-r0;

C = (a1+b1+c1)*t-2*(p1+q1+r1);

D = h*(c0+c1-a1-a0)-V*(p0+q0+r0+p1+q1+r1);

E = h*(c0+c1+a0+a1)-2*V*(r1+r0-p1-p0);

F = h*(c0+c1-a0-a1)-3*V*(r0+r1+p0+p1);

G = (c0+c1+a0+a1)*h-4*V*(r0+r1-p1-p0);

L = t*h*(c1-a1)-h*(r1+r0-p1-p0)-V*t*(p1+q1+r1);

M = t*h*(a1+c1) - h*(r0+r1+p0+p1)-2*V*t*(r1-p1);

H = t^2*h*(c1-a1)-2*t*h*(r1-p1)-V*t^2*(r1+p1+q1);

R =  t*h*(c1-a1)-(r0+r1-p1-p0)*h-3*V*t*(r1+p1);



s  = solve(A, B,C, D, E, F, G, L, M, H, R,   a1, b1, c1,  a0,b0,  c0, r1, r0, q0, q1, p1);

a1 = s.a1
a0 = s.a0
b1 = s.b1
b0 = s.b0
c1 = s.c1
c0 = s.c0
r1 = s.r1
r0 = s.r0
q0 = s.q0
q1 = s.q1
p1 = s.p1




%im = abs(z)
%S30 = simplify(im,'Steps',30)
%S50 = simplify(im,'Steps',50)
%S100 = simplify(im,'Steps',100)
%S500 = simplify(im,'Steps', 500)
%S1000 = simplify(im,'Steps', 1000)

