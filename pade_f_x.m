
syms a1 b1 c1 a0 b0 c0 p0 p1 q0 q1 r0 r1  h t V w z;

p0 = 1;

A = a1+b1+c1+a0+b0+c0;

B = (a1+b1+c1)*t-h*(r1+r0-p1-p0);

C = (a1+b1+c1)*t-2*h*(r1-p1);

D = (c0+c1-a1-a0)-V*(r1+r0-p1-p0);

E = (c0+c1+a0+a1)-V*(r1+r0+p1+p0);

%F = (c0+c1-a0-a1)-V*(r0+r1-p0-p1);

%G = (c0+c1+a0+a1)-V*(r0+r1+p1+p0);

L = t*(c1-a1)-h/2*(r1+r0+p1+p0)-V*t*(r1-p1);

M = t*(a1+c1) - h/3*(r0+r1-p0-p1)-V*t*(r1+p1);

F1 = p1+q1+r1

G1 = p1+q1+r1+p0+q0+r0

H = t^2*h*(c1-a1)-t*h^2*(r1+p1)-V*t^2*h*(r1-p1);

R =  t*(c1-a1)-(r0+r1+p1+p0)*h/4-V*t*(r1-p1);



s  = solve(A, B,C, D, E, F1, G1, L, M, H, R,   a1, b1, c1,  a0,b0,  c0, r1, r0, q0, q1, p1);

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