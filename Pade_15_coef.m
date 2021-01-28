syms a2 b2 c2 a1 b1 c1 a0 b0 c0 p0 p1 q0 q1 b11 r0 r1  h t V w z l;

q0 = 1;


A = a2+b2+c2+a1+b1+c1+a0+b0+c0;
B = p0+q0+r0+p1+q1+r1;
C = (a1+b1+c1)*t+4*t*(a2+b2+c2)-2*h*(r1-p1);

D = (a1+b1+c1)*t+8*t*(a2+b2+c2)-3*h*(r1-p1);

E = (a1+b1+c1)*t+16*t*(a2+b2+c2)-4*h*(r1-p1);


F = (c0+c1+c2-a2-a1-a0)-V*(r0+r1-p0-p1);

G = (c0+c1+c2+a0+a1+a2)-2*V*(r1+r0+p1+p0);

H = (c0+c1+c2-a0-a1-a2)-3*V*(r0+r1-p0-p1);

I = (c0+c1+c2+a0+a1+a2)-4*V*(r0+r1+p1+p0);


K = t*(a1+c1)*h^2+2*t*(a2+c2)*h^2 - h^3/3*(r0+r1-p0-p1)-V*t*h^2*(r1+p1);

L = t*h^3*(c1-a1)+2*t*h^3*(c2-a2)-(r0+r1+p1+p0)*h^4/4-V*t*h^3*(r1-p1);

M = t^2*h*(c1-a1)+4*t^2*h*(c2-a2)-t*h^2*(r1+p1)-V*t^2*h*(r1-p1);

N = t^2*h^2*(c1+a1)+4*t^2*h^2*(c2+a2)-2*t*h^3/3*(r1-p1)-2*V*t^2*h^2*(r1+p1);

O = t^3*h*(c1-a1)+8*t^3*h*(c2-a2)-3*t^2*h^2/2*(r1+p1)-V*t^3*h*(r1-p1);


s  = solve(A,B, C, D, E, F, G,H, I, K, L, M, N, O, a2, b2,c2,  a1, b1, c1,  a0 ,b0, c0, r1, r0,p1,p0, q1);


a2 = s.a2
a1 = s.a1
a0 = s.a0
b2 = s.b2
b1 = s.b1
b0 = s.b0
c2 = s.c2
c1 = s.c1
c0 = s.c0
p1 = s.p1
p0 = s.p0
r1 = s.r1
r0 = s.r0
q1= s.q1


