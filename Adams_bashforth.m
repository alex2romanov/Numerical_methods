n1 = 12;
m1 = 192;
X = linspace(0, pi/2, n1+1);
X = X(1:n1)';
u0 = cos(X)*cos(0);
first = u0;
D = 1;
h = pi/(2*n1);
t = pi/(2*m1);
 for i = 1:m1
     u = zeros(n1, 1);
     u_dop = zeros(n1, 1);
     f = zeros(n1, 1);
     eps = zeros(n1, 1);
     first_con = zeros(n1, 1);
     for j = 1:n1
         
%                 if j == 1
%            
%                     u_dop(j) = (u0(1) + u0(2))/2-t/h*(u0(2)^2/4-u0(1)^2/4); %i шаг по времени, разные по координатам 1 и 0
%                     u_dop0 = (-u0(2)+u0(1))/2-t/h*(u0(1)^2/4-u0(2)^2/4);% j = 0 и -1
%                     u(j) = u0(j) - t/h*(u_dop(j)^2/2-u_dop0^2/2);%j = 0
%             
%                 elseif j == n1
%            
%            
%                     u_dop(j) = (u0(j)+u0(1))/2-t/h*(u0(1)^2/4-u0(j)^2/4);
%                     u(j) = u0(j) - t/h*(u_dop(j)^2/2-u_dop(j-1)^2/2);
%                 else
%              
%                     u_dop(j) = (u0(j)+u0(j+1))/2-t/h*(u0(j+1)^2/4-u0(j)^2/4);
%                     u(j) = u0(j) - t/h*(u_dop(j)^2/2-u_dop(j-1)^2/2);
%                 end
                    
           if j == 1
               u_dop(j) = u0(j)+D*t/(2*h^2)*(u0(2)-2*u0(j)+u0(j+1))+t*u0(j)*(1-u0(j))/2;%для cos(x) на [0;pi/2]
               u_dop0 = u0(2)+D*t/(2*h^2)*(u0(3)-2*u0(2)+u0(j))+t*u0(2)*(1-u0(2))/2;
               u_dop2 = u0(j+1)+D*t/(2*h^2)*(u0(j)-2*u0(j+1)+u0(j+2))+t*u0(j+1)*(1-u0(j+1))/2;
               
               u(j) = u0(j)+D*t/h^2*(u_dop0-2*u_dop(j)+u_dop2)+t*u_dop(j)*(1-u_dop(j));

%                u_dop(j) = u0(j)+D*t/(2*h^2)*(-u0(2)-2*u0(j)+u0(j+1));
%                u_dop0 = -u0(2)+D*t/(2*h^2)*(-u0(3)+2*u0(2)+u0(j));
%                u_dop2 = u0(j+1)+D*t/(2*h^2)*(u0(j)-2*u0(j+1)+u0(j+2));
%                
%                u(j) = u0(j)+D*t/h^2*(u_dop0-2*u_dop(j)+u_dop2);
               
           elseif j == n1
               
               u_dop(j) = u0(j)+D*t/(2*h^2)*(u0(j-1)-2*u0(j)+0)+t*u0(j)*(1-u0(j))/2;
               u_dop0 = u0(j-1)+D*t/(2*h^2)*(u0(j-2)-2*u0(j-1)+u0(j))+t*u0(j-1)*(1-u0(j-1))/2;
               u_dop2 = 0+D*t/(2*h^2)*(u0(j)-2*0-u0(n1))+t*0*(1-0)/2;%только граничное условие для cos [0;pi/2]
               
               u(j) = u0(j)+D*t/h^2*(u_dop0-2*u_dop(j)+u_dop2)+t*u_dop(j)*(1-u_dop(j));

%                u_dop(j) = u0(j)+D*t/(2*h^2)*(u0(j-1)-2*u0(j)+u0(1));
%                u_dop0 = u0(j-1)+D*t/(2*h^2)*(u0(j-2)-2*u0(j-1)+u0(j));
%                u_dop2 = u0(1)+D*t/(2*h^2)*(u0(j)-2*u0(1)+u0(2));
%                
%                u(j) = u0(j)+D*t/h^2*(u_dop0-2*u_dop(j)+u_dop2);
               
   
               
           else
               u_dop(j) = u0(j)+D*t/(2*h^2)*(u0(j-1)-2*u0(j)+u0(j+1))+t*u0(j)*(1-u0(j))/2;
               if j == 2
                   u_dop0 = u0(j-1)+D*t/(2*h^2)*(u0(2)-2*u0(j-1)+u0(j))+t*u0(j-1)*(1-u0(j-1))/2;
                   u_dop2 = u0(j+1)+D*t/(2*h^2)*(u0(j)-2*u0(j+1)+u0(j+2))+t*u0(j+1)*(1-u0(j+1))/2;
               elseif j== n1-1
                   u_dop0 = u0(j-1)+D*t/(2*h^2)*(u0(j-2)-2*u0(j-1)+u0(j))+t*u0(j-1)*(1-u0(j-1))/2;
                   u_dop2 = u0(j+1)+D*t/(2*h^2)*(u0(j)-2*u0(j+1)+0)+t*u0(j+1)*(1-u0(j+1))/2;
               else
               
                    u_dop0 = u0(j-1)+D*t/(2*h^2)*(u0(j-2)-2*u0(j-1)+u0(j))+t*u0(j-1)*(1-u0(j-1))/2;
                    u_dop2 = u0(j+1)+D*t/(2*h^2)*(u0(j)-2*u0(j+1)+u0(j+2))+t*u0(j+1)*(1-u0(j+1))/2;
               end
               
               u(j) = u0(j)+D*t/h^2*(u_dop0-2*u_dop(j)+u_dop2)+t*u_dop(j)*(1-u_dop(j));
               
           end   
        
            
     end
    u_p = zeros(n1+1, 1);
    u_p(1:n1) = u0;
    u_p(n1+1) = 0;
    Y = linspace(0, pi/2, n1+1);
    if i <=192
        hold on
        plot(Y, u_p);
  
    end
    u0 = u;
    
   
    
 end
 hold off
        