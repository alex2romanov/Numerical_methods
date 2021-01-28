syms a2 b2 c2 a1 b1 c1 a0 b0 c0 p0 p1 q0 q1 r0 r1 r2 p2 q2  h t V w z l;


q2 = 0;
p2 = -1;
r2 = 1;

A = a2+b2+c2+a1+b1+c1+a0+b0+c0;

B = (a1+b1+c1)*t+2*t*(a2+b2+c2)-p0-q0-r0-p1-q1-r1-p2-q2-r2;

C = (a1+b1+c1)*t+4*t*(a2+b2+c2)-2*(p1+q1+r1)-4*(p2+q2+r2);

D = (a1+b1+c1)*t+8*t*(a2+b2+c2)-3*(r1+p1+q1)-12*(r2+p2+q2);

E = (a1+b1+c1)*t+16*t*(a2+b2+c2)-4*(r1+p1+q1)-32*(r2+p2+q2);


F = (c0+c1+c2-a2-a1-a0)*h-V*(p0+q0+r0+p1+q1+r1+p2+q2+r2);

G = (c0+c1+c2+a0+a1+a2)*h-2*V*(r2+r1+r0-p2-p1-p0);

H = (c0+c1+c2-a0-a1-a2)*h-3*V*(r0+r1+r2+p2+p1+p0);

I = (c0+c1+c2+a0+a1+a2)*h-4*V*(r0+r1+r2-p2-p1-p0);

J = (c1-a1)*t*h+2*t*h*(c2-a2)-h*(r0+r1+r2-p0-p1-p2)-V*t*(p1+q1+r1)-2*V*t*(p2+q2+r2);

K = t*(a1+c1)*h^2+2*t*(a2+c2)*h^2 - h^2*(r0+r1+r2+p0+p1+p2)-2*V*t*h*(r1-p1)-4*V*t*h*(r2-p2);

L = t*h^3*(c1-a1)+2*t*h^3*(c2-a2)-(r0+r1+r2-p2-p1-p0)*h^3-3*V*t*h^2*(r1+p1)-6*V*t*h^2*(r2+p2);

M = t^2*h*(c1-a1)+4*t^2*h*(c2-a2)-2*t*h*(r1-p1)-4*t*h*(r2-p2)-V*t^2*(r1+p1+q1)-4*V*t^2*(r2+p2+q2);


N = t^2*h^2*(c1+a1)+4*t^2*h^2*(c2+a2)-2*t*h^2*(r1+p1)-4*t*h^2*(r2+p2)-2*V*t^2*h*(r1-p1)-8*V*t^2*h*(r2-p2);

O = t^3*h*(c1-a1)+8*t^3*h*(c2-a2)-3*t^2*h*(r1-p1)-12*t^2*h*(r2-p2)-V*t^3*(r1+p1+q1)-8*V*t^3*(r2+p2+q2);

%P = (c0+c1+c2-a0-a1-a2)*h-5*V*(r0+r1+r2+p0+p1+p2);

%R = (c1+a1)*t*h^4+2*t*h^4*(c2+a2)-h^4*(r0+r1+r2+p0+p1+p2)-4*V*t*h^3*(r1-p1)-8*V*t*h^3*(r2-p2);

%RR = t^2*h^3*(c1-a1)+4*t^2*h^3*(c2-a2)-2*t*h^3*(r1-p1)-4*t*h^3*(r2-p2)-3*V*t^2*h^2*(r1+p1)-12*V*t^2*h^2*(r2+p2);

s  = solve(A,B, C, D, E, F, G,H, I,J,  K, L, M, N, O, a2, b2,c2, a0,b0, c0,  a1, b1, c1, r1, r0,p1,p0, q1, q0);


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
q1 = s.q1
q0 = s.q0



