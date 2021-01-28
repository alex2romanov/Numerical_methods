format long
n1 = 48;
m1 = 192;
V = 1;
% надо учесть что Курант tay*max|u|/h
X = linspace(0, 2*pi, n1+1);
Y = X
X = X(1:n1)';
u0 = sin(X)*cos(0);
first = u0

%u0 = ones(n1, 1);
%[result, abc] = solve_diff(n1, m1, V, u0);
%for i = 1:m1+1
    
%    plot(X, result(:, i))
%    hold on
%end
U = zeros(n1, m1+1);
F = zeros(n1, 2);
epsilon = zeros(n1, m1+1);


U(:, 1) = u0;
N = [10, 20, 40, 80, 160, 320, 640];
M = [20, 40, 80, 160, 320, 640];
%N = [5, 24];
%M = [5, 24];
norm = zeros(length(N), length(M));

for l = 1:length(N)
    for k = 1:length(M)
        
        X = linspace(0, 2*pi, N(l)+1);
        X = X(1:N(l))';
        u0 = sin(X)*cos(0);
        first = u0;
        h = 2*pi/N(l);
        t = 2*pi/M(k);
        
        for i = 1:M(k)
            u = zeros(N(l), 1);
            u_dop = zeros(N(l), 1);
            f = zeros(N(l), 1);
            eps = zeros(N(l), 1);
            first_con = zeros(N(l), 1);
            
            for j = 1:N(l)
        
                if j == 1
                    u_dop(j) = (u0(1) + u0(2))/2-t/h*(u0(2)^2/4-u0(1)^2/4-cos(h)*sin(i*t)/2+cos(2*h)*cos(i*t)^2/8+cos(0)*sin(i*t)/2-cos(2*0)*cos(i*t)^2/8); %i шаг по времени, разные по координатам 1 и 0
                    u_dop0 = (-u0(2)+u0(1))/2-t/h*(u0(1)^2/4-u0(2)^2/4-cos(0*h)*sin(i*t)/2+cos(0*h)*cos(i*t)^2/8+cos(-1*h)*sin(i*t)/2-cos(-2*h)*cos(i*t)^2/8);% j = 0 и -1
                    u(j) = u0(j) - t/h*(u_dop(j)^2/2-u_dop0^2/2-cos((1/2)*h)*sin((i+1/2)*t)+cos((1)*h)*cos((i+1/2)*t)^2/4  +  cos((-1/2)*h)*sin((i+1/2)*t)-cos((-1)*h)*cos((i+1/2)*t)^2/4);%j = 0
           
                    %u_dop(j) = (u0(1) + u0(2))/2-t/h*(u0(2)^2/4-u0(1)^2/4); %i шаг по времени, разные по координатам 1 и 0
                    %u_dop0 = (-u0(2)+u0(1))/2-t/h*(u0(1)^2/4-u0(2)^2/4);% j = 0 и -1
                    %u(j) = u0(j) - t/h*(u_dop(j)^2/2-u_dop0^2/2);%j = 0
            
                elseif j == N(l)
            
                    u_dop(j) = (u0(j)+u0(1))/2-t/h*(u0(1)^2/4-u0(j)^2/4-cos(j*h)*sin(i*t)/2+cos(2*j*h)*cos(i*t)^2/8+cos((j-1)*h)*sin(i*t)/2-cos(2*(j-1)*h)*cos(i*t)^2/8);
                    u(j) = u0(j) - t/h*(u_dop(j)^2/2-u_dop(j-1)^2/2-cos((j+1/2-1)*h)*sin((i+1/2)*t)+cos((2*j-1)*h)*cos((i+1/2)*t)^2/4+cos((j-1/2-1)*h)*sin((i+1/2)*t)-cos((2*j-3)*h)*cos((i+1/2)*t)^2/4);
           
                    %u_dop(j) = (u0(j)+u0(1))/2-t/h*(u0(1)^2/4-u0(j)^2/4);
                    %u(j) = u0(j) - t/h*(u_dop(j)^2/2-u_dop(j-1)^2/2);
                else
                    u_dop(j) = (u0(j)+u0(j+1))/2-t/h*(u0(j+1)^2/4-u0(j)^2/4-cos(j*h)*sin(i*t)/2+cos(2*j*h)*cos(i*t)^2/8+cos((j-1)*h)*sin(i*t)/2-cos(2*(j-1)*h)*cos(i*t)^2/8);
                    u(j) = u0(j) - t/h*(u_dop(j)^2/2-u_dop(j-1)^2/2-cos((j+1/2-1)*h)*sin((i+1/2)*t)+cos((2*j-1)*h)*cos((i+1/2)*t)^2/4+cos((j-1/2-1)*h)*sin((i+1/2)*t)-cos((2*j-3)*h)*cos((i+1/2)*t)^2/4);
             
                    %u_dop(j) = (u0(j)+u0(j+1))/2-t/h*(u0(j+1)^2/4-u0(j)^2/4);
                    %u(j) = u0(j) - t/h*(u_dop(j)^2/2-u_dop(j-1)^2/2);
                end
        
            
            end

            %for j = 1:n1
            %        first_con(j) = sin((j-1)*h)*cos(i*t);
            %end
    
            u_p = zeros(N(l)+1, 1);
            u_p(1:N(l)) = u0;
            u_p(N(l)+1) = 0;
            
            %norm(1, i)= max(abs(u0 - first));
    
    
            %if ((N(l) == 24) && (M(k) == 48))
            %    if i <= 10
            %        hold on
            %        plot(Y, u_p);
 
            %    end
            %end
            
            
%             xlabel('h координата')
    
             % first iteration with epsilons

            for j = 1:N(l)
                if j == 1
                    f(j) = -u0(2)+4*u0(j)+u0(j+1)+u(2)-4*u(j)-u(j+1)+3*t/(4*h)*(u0(2)^2+u(2)^2-u0(j+1)^2-u(j+1)^2)+3*t/(2*h)*(-cos(-h)*sin(i*t)+cos(-2*h)*cos(i*t)^2/4+cos(h)*sin(i*t)-cos(2*h)*cos(i*t)^2/4-cos(-1*h)*sin((i+1)*t)+cos(-2*h)*cos((i+1)*t)^2/4+cos(1*h)*sin((i+1)*t)-cos(2*h)*cos((i+1)*t)^2/4);
                    %+3*t/(2*h)*(-cos(-h)*sin(i*t)+cos(-2*h)*cos(i*t)^2/4+cos(h)*sin(i*t)-cos(2*h)*cos(i*t)^2/4-cos(-1*h)*sin((i+1)*t)+cos(-2*h)*cos((i+1)*t)^2/4+cos(1*h)*sin((i+1)*t)-cos(2*h)*cos((i+1)*t)^2/4);% вычитание единицы из-за специфики матлаба , -u(2) работает только если u(1)=0. 
                elseif j == N(l)
                    f(j) = u0(j-1)+4*u0(j)+u0(1)-u(j-1)-4*u(j)-u(1)+3*t/(4*h)*(u0(j-1)^2+u(j-1)^2-u0(1)^2-u(1)^2)+3*t/(2*h)*(-cos((j-2)*h)*sin(i*t)+cos(2*(j-2)*h)*cos(i*t)^2/4+cos((j)*h)*sin(i*t)-cos(2*(j)*h)*cos(i*t)^2/4-cos((j-2)*h)*sin((i+1)*t)+cos(2*(j-2)*h)*cos((i+1)*t)^2/4+cos((j)*h)*sin((i+1)*t)-cos(2*(j)*h)*cos((i+1)*t)^2/4);% также пользуюсь, что u(1)=0
                else
                    f(j) = u0(j-1)+4*u0(j)+u0(j+1)-u(j-1)-4*u(j)-u(j+1)+3*t/(4*h)*(u0(j-1)^2+u(j-1)^2-u0(j+1)^2-u(j+1)^2)+3*t/(2*h)*(-cos((j-2)*h)*sin(i*t)+cos(2*(j-2)*h)*cos(i*t)^2/4+cos((j)*h)*sin(i*t)-cos(2*(j)*h)*cos(i*t)^2/4-cos((j-2)*h)*sin((i+1)*t)+cos(2*(j-2)*h)*cos((i+1)*t)^2/4+cos((j)*h)*sin((i+1)*t)-cos(2*(j)*h)*cos((i+1)*t)^2/4);
                    %+3*t/(2*h)*(-cos((j-2)*h)*sin(i*t)+cos(2*(j-2)*h)*cos(i*t)^2/4+cos((j)*h)*sin(i*t)-cos(2*(j)*h)*cos(i*t)^2/4-cos((j-2)*h)*sin((i+1)*t)+cos(2*(j-2)*h)*cos((i+1)*t)^2/4+cos((j)*h)*sin((i+1)*t)-cos(2*(j)*h)*cos((i+1)*t)^2/4)
                end
            end

            E = diag(ones(N(l), 1)*4);

            for j = 1:N(l)-1

                if j == 1
                    E(2, 1) = 1-3*t/(2*h)*u(j);
                else
                    E(j+1, j)= 1-3*t/(2*h)*u(j);% крайнее значение u(n1-2) и +1 сдвиг
                end

            end

            for j = 2:N(l)

                if j == N(l)
                    E(j-1, j) = 1+3*t/(2*h)*u(j);
                else
                    E(j-1, j) = 1+3*t/(2*h)*u(j);
                end

            end

            E(1, N(l))=1-3*t/(2*h)*u(N(l));

            E(N(l), 1)=1+3*t/(2*h)*u(1);


            eps = inv(E)*f;



            u0 = u+eps;
            %u0 = u;
            %norm(1, i)= max(abs(u - first_con));
    
%     
%     
%     
%     %  second iteration with epsilons
    
%     for j = 1:n1
%         if j == 1
%             f(j) = -u0(2)+4*u0(j)+u0(j+1)+u(2)-4*u(j)-u(j+1)+3*t/(4*h)*(u0(2)^2+u(2)^2-u0(j+1)^2-u(j+1)^2)+3*t/(2*h)*(-cos(-h)*sin(i*t)+cos(-2*h)*cos(i*t)^2/4+cos(h)*sin(i*t)-cos(2*h)*cos(i*t)^2/4-cos(-1*h)*sin((i+1)*t)+cos(-2*h)*cos((i+1)*t)^2/4+cos(1*h)*sin((i+1)*t)-cos(2*h)*cos((i+1)*t)^2/4);% вычитание единицы из-за специфики матлаба , -u(2) работает только если u(1)=0. 
%         elseif j == n1
%             f(j) = u0(j-1)+4*u0(j)+u0(1)-u(j-1)-4*u(j)-u(1)+3*t/(4*h)*(u0(j-1)^2+u(j-1)^2-u0(1)^2-u(1)^2)+3*t/(2*h)*(-cos((j-2)*h)*sin(i*t)+cos(2*(j-2)*h)*cos(i*t)^2/4+cos((j)*h)*sin(i*t)-cos(2*(j)*h)*cos(i*t)^2/4-cos((j-2)*h)*sin((i+1)*t)+cos(2*(j-2)*h)*cos((i+1)*t)^2/4+cos((j)*h)*sin((i+1)*t)-cos(2*(j)*h)*cos((i+1)*t)^2/4);% также пользуюсь, что u(1)=0
%         else
%             f(j) = u0(j-1)+4*u0(j)+u0(j+1)-u(j-1)-4*u(j)-u(j+1)+3*t/(4*h)*(u0(j-1)^2+u(j-1)^2-u0(j+1)^2-u(j+1)^2)+3*t/(2*h)*(-cos((j-2)*h)*sin(i*t)+cos(2*(j-2)*h)*cos(i*t)^2/4+cos((j)*h)*sin(i*t)-cos(2*(j)*h)*cos(i*t)^2/4-cos((j-2)*h)*sin((i+1)*t)+cos(2*(j-2)*h)*cos((i+1)*t)^2/4+cos((j)*h)*sin((i+1)*t)-cos(2*(j)*h)*cos((i+1)*t)^2/4);
%         end
%     end
%     
%     E = diag(ones(n1, 1)*4);
%     
%     for j = 1:n1-1
%         
%         if j == 1
%             E(2, 1) = 1-3*t/(2*h)*u(j);
%         else
%             E(j+1, j)= 1-3*t/(2*h)*u(j);% крайнее значение u(n1-2) и +1 сдвиг
%         end
%         
%     end
%    
%     for j = 2:n1
%         
%         if j == n1
%             E(j-1, j) = 1+3*t/(2*h)*u(j);
%         else
%             E(j-1, j) = 1+3*t/(2*h)*u(j);
%         end
%         
%     end
%         
%     E(1, n1)=1-3*t/(2*h)*u(n1);
%     
%     
%     E(n1, 1)=1+3*t/(2*h)*u(1);
%     
%     
%     
%   
%     
%     eps = inv(E)*f;
%     
%     if i <= 1
%         hold on
%         plot(X, eps);
%         
%     end
%     
%     u0 = u+eps;
%     if i == 1
%         F(:, i+1) = f;
%     end
%     norm(2, i)= max(abs(u0 - first_con));
%     
%     
           

%             if ((N(l) == 192) && (M(k) == 384))
%                 if i <=200
%                     hold on
%                     plot(X, u0);
%  
%                 end
%             end


        end
        norm(l, k) = log(max(abs(first - u0)));
    end
end
hold off

%otv = max(abs(first - u0)) % тут записано значение нормы
%plot(X(1:10), norm(3, 1:10));

contour(M, N, norm, 'ShowText', 'on')
xlabel('M')
ylabel('N')