function U = diff_krank(N, M, D, u0, nu)

X = linspace(0, 2*pi, N+1);


% a0 = 2*(6*D*t/(h*h)-1);
% a1 = 2*(6*D*t/(h*h)+1);
% b0 = -4*(6*D*t/(h*h)+5);
% b1 = -4*(6*D*t/(h*h)-5);
% p0 = -t;
% p1 = -t;
% q0 = -10*t;
% q1 = -10*t;

a0 = 2*(6*nu-1);
a1 = 2*(6*nu+1);
b0 = -4*(6*nu+5);
b1 = -4*(6*nu-5);


% a1 = -V*t/h; % Коэффициенты для Кранк- Николсона
% b1 = 4;
% c1 = V*t/h;
% 
% a0 = -V*t/h;
% b0 = -4;
% c0 = V*t/h;
% q0 = 0;
% q1 = 0;% правая часть df/dx
% p0 = -t/(4*h);  
% r0 = t/(4*h);
% p1 = -t/(4*h);
% r1 = t/(4*h);



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








A = diag(ones(N, 1)*b0)+diag(ones(N-1,1)*a0,1)+diag(ones(N-1,1)*a0,-1);
A(1, N)=a0; A(N, 1)=a0;




A0 = diag(ones(N, 1)*(-b1))+diag(ones(N-1,1)*(-a1),1)+diag(ones(N-1,1)*(-a1),-1);
A0(1, N)=-a1; A0(N, 1)=-a1;


% C0 = diag(ones(N, 1)*q0)+diag(ones(N-1,1)*r0,1)+diag(ones(N-1,1)*p0,-1);
% C1 = diag(ones(N, 1)*q1)+diag(ones(N-1,1)*r1,1)+diag(ones(N-1,1)*p1,-1);
% 
% 
% C0(1, N)=p0; C0(N, 1)=r0;
% C1(1, N)=p1; C1(N, 1)=r1;

%C0 = sparse(C0);
%C1 = sparse(C1);


Ainv = inv(A);


u = u0;
%U = zeros(N+1, M+1);
U(:, 1) = u0;


% f1 = zeros(N, 1);
% f0 = zeros(N, 1);



for j = 1:M
    u_plot = zeros(N+1);
    u_plot(1:N) = u;
    cur = A0*u;
    plot(X, u_plot)
    u = A\cur;
    U(:, j+1) = u;
    
    hold on
              
            
end 
hold off



% for j = 1:M % тут считается шаг по времени
%             
%     for k = 1:N % чтобы посчитать начальные значения по f
%         
%           f1(k) = cos((k-1)*h)*sin(t*j)+V*sin((k-1)*h)*cos(t*j);
%           f0(k) = cos((k-1)*h)*sin(t*(j-1))+V*sin((k-1)*h)*cos(t*(j-1));
% 
%           
%     end
% 
%     cur = A0*u +C0*f0+C1*f1; % Задание реккурентной формулы
%     u = Ainv*cur; 
%     U(:, j+1) = u;
% 
%               
% end


end 

% Итого точность 2