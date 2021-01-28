function U = solve_diff(N, M, V, u0)

h = 2*pi/N;
t = 2*pi/M;

%a1 = -V*t/h;
%b1 = 4;
%c1 = V*t/h;

%a0 = -V*t/h;
%b0 = -4;
%c0 = V*t/h;



a0 = -(2*h+3*V*t)/(h*t);
b0 = -8/t;
c0 = -(2*h-3*V*t)/(h*t);

a1 = (2*h-3*V*t)/(h*t);
b1 = 8/t;
c1 = (2*h+3*V*t)/(h*t);




%a1 = (2*h-3*V*t)/(h*t);
%c1 = (2*h+3*V*t)/(h*t);
%a0 = (-2*h-3*V*t)/(h*t);
%b0 = -8/t;
%c0 = (-2*h+3*V*t)/(h*t);
%b1=8/t; 
%p0 = 1;
%q0 = 4;
%r0 = 1;
%p1 = 1;
%q1 = 4;
%r1 = 1;




A = diag(ones(N, 1)*b1)+diag(ones(N-1,1)*c1,1)+diag(ones(N-1,1)*a1,-1);
A(1, N)=a1; A(N, 1)=c1;

A = sparse(A);
B = diag(ones(N, 1)*(-b0))+diag(ones(N-1,1)*(-c0),1)+diag(ones(N-1,1)*(-a0),-1);
B(1, N)=-a0; B(N, 1)=-c0;

%C0 = diag(ones(N, 1)*q0)+diag(ones(N-1,1)*p0,1)+diag(ones(N-1,1)*r0,-1);
%C1 = diag(ones(N, 1)*q1)+diag(ones(N-1,1)*p1,1)+diag(ones(N-1,1)*r1,-1);
%C0(1, N)=p0; C0(N, 1)=r0;
%C1(1, N)=p1; C1(N, 1)=r1;
%C0 = sparse(C0);
%C1 = sparse(C1);

B = sparse(B);
Ainv = inv(A);





%X = linspace(0, 2*pi, N+1);
%X = X(1:N)';

u = u0;
U = zeros(N, M+1);
U(:, 1) = u0;
f1 = zeros(N, 1);
f0 = zeros(N, 1);





for j = 1:M
            
    cur = B*u;
    u = A\cur;
    
    U(:, j+1) = u;
              
            
end 

end 