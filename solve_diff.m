function [res, abc] = solve_diff(N, M, V, u0)

h = 2*pi/N;
t = 2*pi/M;
mu = V*t/h;

%a1 = -V*t/h; % Коэффициенты для Кранк- Николсона
%b1 = 4;
%c1 = V*t/h;

%a0 = -V*t/h;
%b0 = -4;
%c0 = V*t/h;



%a0 = -(2*h+3*V*t)/(h*t);
%b0 = -8/t;
%c0 = -(2*h-3*V*t)/(h*t); % Какие-то из этих двух коэффициенты для компакта

%a1 = (2*h-3*V*t)/(h*t);
%b1 = 8/t;
%c1 = (2*h+3*V*t)/(h*t);



%a1 = (2*h-3*V*t)/(h*t);
%c1 = (2*h+3*V*t)/(h*t);
%a0 = (-2*h-3*V*t)/(h*t);
%b0 = -8/t;
%c0 = (-2*h+3*V*t)/(h*t);
%b1=8/t; 



%q0 = t/2;
%q1 = t/2
%p0 = 0;
%r0 = 0; q0 = t/2 и q1 = t/2 коэффициенты для правой части f для КН
%p1= 0;
%r1 = 0;




%q0 = 0
%q1 = 0
%p0 = -t/(4*h);  Коэффициенты для правой части df/dx
%r0 = t/(4*h);
%p1 = -t/(4*h);
%r1 = t/(4*h);




%a0 = (-4*V^2*t^2-13*V*h*t-9*h^2)/(2*h*t);
%a1 = (2*V^2*t^2-8*V*h*t+4*h^2)/(h*t);
%b0 = (4*V^2*t^2-18*h^2)/(h*t);
%b1 = (8*h-4*V*t)*(2*h+V*t)/(h*t);
%c0 = (-4*V^2*t^2+13*V*h*t-9*h^2)/(2*h*t);
%c1 = (2*V^2*t^2+8*V*h*t+4*h^2)/(h*t);
%r0=(2*h-V*t);
%r1 = (3*h+V*t);
%p0 = 2*h+V*t;
%p1 = (3*h-V*t);
%q0 = 10*h;
%q1 = 11*h;
%q11 = -h;
%a2 = (h-V*t)/(2*t);
%b2 = (2*h)/t;
%c2 = (h+V*t)/(2*t);


a0 = -2+2*mu;
a1 = 4;
a2 = -2+2*mu;

b0 = 4*mu;
b1 = 0;
b2 = -4*mu;

c0 = 2-2*mu;
c1 = -4;
c2 = 2+2*mu;

r0 = -t;
r1 = 0;
r2 = t;
p0 = t;
p1 = 0;
p2 = -t;
q0 = 0;
q1 = 0;
q2 = 0;


A = diag(ones(N, 1)*b1)+diag(ones(N-1,1)*c1,1)+diag(ones(N-1,1)*a1,-1);
A(1, N)=a1; A(N, 1)=c1;

A2 = diag(ones(N, 1)*b2)+diag(ones(N-1,1)*c2,1)+diag(ones(N-1,1)*a2,-1);
A2(1, N) = a2; A2(N, 1) = c2;
A2(1,:)= ones(N, 1);



A = sparse(A);

A0 = diag(ones(N, 1)*(b0))+diag(ones(N-1,1)*(c0),1)+diag(ones(N-1,1)*(a0),-1);
A0(1, N)=a0; A0(N, 1)=c0;

A1 = diag(ones(N, 1)*(b1))+diag(ones(N-1,1)*(c1),1)+diag(ones(N-1,1)*(a1),-1);
A1(1, N)=a1; A1(N, 1)=c1;% потом не нужно оборачивать A1
A1 = sparse(A1);

C0 = diag(ones(N, 1)*q0)+diag(ones(N-1,1)*r0,1)+diag(ones(N-1,1)*p0,-1);
C1 = diag(ones(N, 1)*q1)+diag(ones(N-1,1)*r1,1)+diag(ones(N-1,1)*p1,-1);

%A11 = diag(ones(N, 1)*a11);

C0(1, N)=p0; C0(N, 1)=r0;
C1(1, N)=p1; C1(N, 1)=r1;

C0 = sparse(C0);
C1 = sparse(C1);

%A11 = sparse(A11);

%A0 = sparse(A0);
Ainv = inv(A);
Ainv_2 = inv(A2);


%w = zeros(48, 48);
%w(1:24, 1:24) = Ainv_2;
%w(1:24, 25:48) = Ainv_2;
%w(25:48, 1:24) = eye(24);
%w(25:48, 25:48) = zeros(24, 24);



%w(1:16, 1:16)= -inv(A2)*A1;
%w(1:16, 17:32) = -inv(A2)*A0;
%w(1:16, 33:48) = -inv(A2)*A11;
%w(17:32,1:16) = eye(16);
%w(17:32, 17:32) = zeros(16, 16);
%w(17:32, 33:48) = zeros(16, 16);
%w(33:48, 1:16) = zeros(16, 16);
%w(33:48, 17:32) = eye(16);
%w(33:48, 33:48) = zeros(16, 16);

%res = eig(w);
abc = 2;
res = -Ainv_2*A1;



%X = linspace(0, 2*pi, N+1);
%X = X(1:N)';

u = u0;
U = zeros(N, M+1);
P = zeros(5, M+1);
U(:, 1) = u0;
f1 = zeros(N, 1);
f0 = zeros(N, 1);
f11 = zeros(N, 1);
L = zeros(1, M+1);
ftest = zeros(N, 1);



for j = 1:M % тут считается шаг по времени
     for k = 1:N % чтобы посчитать начальные значения по f
        
          %f1(k) = V*sin((k-1)*h)*cos(t*j)+cos((k-1)*h)*sin(t*j);
          %f0(k) = V*sin((k-1)*h)*cos(t*(j-1))+cos((k-1)*h)*sin(t*(j-1));
          %f11(k) = V*sin((k-1)*h)*cos(t*(j-2))+cos((k-1)*h)*sin(t*(j-2));
          
          
          %f1(k) = V*sin((k-1)*h)*cos(t*j);
          %f0(k) = V*sin((k-1)*h)*cos(t*(j-1));
          %f11(k) = V*sin((k-1)*h)*cos(t*(j-2));
          
    end 

    if j == 1
        %cur = B*u +C0*f0+C1*f1;
        %l = Ainv*cur; % u^1
        
        %U(:, j+1) = l;
        
        for k = 1:N 
            
          %ftest(k) = V*sin((k-1)*h)*cos(t*j);
          ftest(k) = 1;

        end
        
        U(:, j+1) = ftest;
    
    else
        
        cur1 = A1*ftest + A0*u + C0*f0+C1*f1+C11*f11; % то есть сначала тут ftest = u1, u = u0
        
        %first = A1*ftest;
        %second = A0*u;
        %third = C0*f0;
        %four = C1*f1;
        %five = C11*f11;
        
        %P(1, j+1) = first(1);
        %P(2, j+1) = second(1);
        %P(3, j+1) = third(1);
        %P(4, j+1) = four(1);
        %P(5, j+1) = five(1);
        
        u = ftest; %u = u1
        ftest = Ainv_2\cur1; %ftest = u2
        U(:, j+1) = ftest;
        %L(j+1) = first(1);
        
    end
    
    
              
end


end 