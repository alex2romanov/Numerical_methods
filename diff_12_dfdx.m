syms a0 b0 c0 a1 b1 c1 p0 p1 q0 q1 r0 r1 D h t;

A0 = q0 - 1;


A = a0+b0+c0+a1+b1+c1;

B = p0+q0+r0+p1+q1+r1;

V = p0+q0+r0;

C = (a0+b0+c0)*t-(r0+r1-p0-p1)*h;


Q = (a0+b0+c0)*t^2-2*h*t*(r0-p0);

E = (a0+b0+c0)*t^3-3*h*t^2*(r0-p0);

F = (c0+c1-a0-a1);

G = (c0+c1+a0+a1)*h^2+2*D*h*(r0+r1-p0-p1);

H = (c0+c1-a0-a1)*h^3+3*D*h^2*(r0+r1+p0+p1);

L = t*h*(c0-a0)-(r0+r1+p0+p1)*h^2/2;

M = (c0+a0)*h^2*t-(r0+r1-p0-p1)*h^3/3+2*D*h*t*(r0-p0);

%S = (c0-a0)*h*t^2-(r0+p0)*t*h^2;

P = (c0+c1+a0+a1)*h^4+4*D*(r0+r1-p0-p1)*h^3;


s  = solve(A0, A,B, V,C, Q,E,  F, G,H, L, M,P, a1, b1, c1,a0, b0,c0, p1,q1, r1,p0,  q0, r0);


a1 = s.a1
b1 = s.b1
c1 = s.c1
a0 = s.a0
b0 = s.b0
c0 = s.c0
p1 = s.p1
q1 = s.q1
r1 = s.r1
p0 = s.p0
q0 = s.q0
r0 = s.r0
