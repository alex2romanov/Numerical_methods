N1 = [10];
M1 = [10];
V = 1;
Z = zeros(length(N1), length(M1));
T = zeros(length(N1), length(M1));
mass=0;
pech = 0;
pech1 = 0;
upr = 0;
for nn = 1:length(N1)
    for mm = 1:length(M1)
        
        tic % timer
        N = N1(nn);
        M = M1(mm);
        h = 2*pi/N;
        t = 2*pi/M;
        
        %a1 = 2-3*V*t/h;
        %c1 = 2+3*V*t/h;
        %a0 = -2-3*V*t/h;
        %b0 = -8;
        %c0 = -2+3*V*t/h;
        %b1=8; 
        %p0 = 1;
        %q0 = 0;
        %r0 = -1;
        %p1 = 1;
        %q1 = 0;
        %r1 = -1;
        
        a1 = (2*h-3*V*t)/(h*t);
        c1 = (2*h+3*V*t)/(h*t);
        a0 = (-2*h-3*V*t)/(h*t);
        b0 = -8/t;
        c0 = (-2*h+3*V*t)/(h*t);
        b1=8/t; 
        p0 = 1;
        q0 = 4;
        r0 = 1;
        p1 = 1;
        q1 = 4;
        r1 = 1;
        

        A = diag(ones(N, 1)*b1)+diag(ones(N-1,1)*c1,1)+diag(ones(N-1,1)*a1,-1);
        A(1, N)=a1; A(N, 1)=c1;

        B = diag(ones(N, 1)*(-b0))+diag(ones(N-1,1)*(-c0),1)+diag(ones(N-1,1)*(-a0),-1);
        B(1, N)=-a0; B(N, 1)=-c0;
        C0 = diag(ones(N, 1)*q0)+diag(ones(N-1,1)*p0,1)+diag(ones(N-1,1)*r0,-1);
        C1 = diag(ones(N, 1)*q1)+diag(ones(N-1,1)*p1,1)+diag(ones(N-1,1)*r1,-1);
        C0(1, N)=p0; C0(N, 1)=r0;
        C1(1, N)=p1; C1(N, 1)=r1;
        
        %B = sparse(B);
        Ainv = inv(A);
        
        u0 = zeros(N, 1);
        u_etal = zeros(N, 1);
        
        for i = 1:N
            u0(i)=sin((i-1)*h);
            u_etal(i) = sin((i-1)*h);
        end
       
        
        f1 = zeros(N, 1);
        f0 = zeros(N, 1);
        for j = 1:M % тут считается шаг по времени
            
            for k = 1:N % чтобы посчитать начальные значения по f
                f1(k) = V*cos((k-1)*h)*cos(t*j)-sin((k-1)*h)*sin(t*j);
                f0(k) = V*cos((k-1)*h)*cos(t*(j-1))-sin((k-1)*h)*sin(t*(j-1));
            end
             
            cur = B*u0 +C0*f0+C1*f1;
            u0 = Ainv*cur;   
              
        end
        
        %upr(length(upr)+1) = u0(1);
        
        %x = linspace(1, N+1, N+1);
        %if (N == 10 && M == 10)
        %    plot(x, upr)
        %end

        
        %if (N == 10 && M == 20)
        %    plot(x-1, u0)
        %    hold off
        %end
       
    
        %norma1 = sqrt(sum((u0-u_etal).*(u0-u_etal)).*h);
        %Z(nn, mm) = norma1;
        %T(nn, mm) = toc;
        
        x = linspace(1, N+1, N+1);
        if (M == 20 && N == 10)
            arr = zeros(1, length(N)+1);
            arr(1) = (u0(length(u0)));
            for i = 1:length(u0)
                arr(i+1) = u0(i);
            end
            plot(x-1, arr)
            hold on
        end
        
        if (M == 40 && N == 10)
            arr = zeros(1, length(N)+1);
            arr(1) = (u0(length(u0)));
            for i = 1:length(u0)
                arr(i+1) = u0(i);
            end
            plot(x-1, arr)
            hold on
            
        end
        
        if (M == 80 && N == 10)
            arr = zeros(1, length(N)+1);
            arr(1) = (u0(length(u0)));
            for i = 1:length(u0)
                arr(i+1) = u0(i);
            end
            plot(x-1, arr)
            hold on
            
        end
        
        if (M == 80 && N == 20)
            
            arr = zeros(1, length(N)+1);
            arr(1) = (u0(length(u0)));
            for i = 1:length(u0)
                arr(i+1) = u0(i);
            end
            
            %arr1 = zeros(1, 11);
            %for i = 1:5
            %    arr1(i) = (arr(2*i)+arr(2*i-1))/2;
            %end
            %arr1(i)= arr(11);
            %for i = 6:10
            %    arr1(i) = (arr(2*i)+arr(2*i+1))/2;
            %end
            x = linspace(1, 11, 11);
            %plot(x-1, arr1);
            %hold on
            
            n=20; % количество меток
            step=ceil(n/10);
            disp(length(arr(1:step:end)));
            plot(x-1,arr(1:step:end),'--')

            hold on

        end
        
        if (M == 80 && N == 20)
            
            arr = zeros(1, length(N)+1);
            arr(1) = (u0(length(u0)));
            for i = 1:length(u0)
                arr(i+1) = u0(i);
            end
            
            %arr1 = zeros(1, 11);
            %for i = 1:5
            %    arr1(i) = (arr(2*i)+arr(2*i-1))/2;
            %end
            %arr1(i)= arr(11);
            %for i = 6:10
            %    arr1(i) = (arr(2*i)+arr(2*i+1))/2;
            %end
            x = linspace(1, 11, 11);
            %plot(x-1, arr1);
            %hold on
            
            n=20; % количество меток
            step=ceil(n/10);
            disp(length(arr(1:step:end)));
            plot(x-1,arr(1:step:end),'--')

            hold on

         end
        
         if (M == 80 && N == 40)
            
            arr = zeros(1, length(N)+1);
            arr(1) = (u0(length(u0)));
            for i = 1:length(u0)
                arr(i+1) = u0(i);
            end
            
            %arr1 = zeros(1, 11);
            %for i = 1:5
            %    arr1(i) = (arr(2*i)+arr(2*i-1))/2;
            %end
            %arr1(i)= arr(11);
            %for i = 6:10
            %    arr1(i) = (arr(2*i)+arr(2*i+1))/2;
            %end
            x = linspace(1, 11, 11);
            %plot(x-1, arr1);
            %hold on
            
            n=40; % количество меток
            step=ceil(n/10);
            plot(x-1,arr(1:step:end),'*')

            hold on

        end
            
    end
end

for j = 1:M % тут считается шаг по времени
            
    for k = 1:N % чтобы посчитать начальные значения по f
          f1(k) = V*cos((k-1)*h)*cos(t*j)-sin((k-1)*h)*sin(t*j);
          f0(k) = V*cos((k-1)*h)*cos(t*(j-1))-sin((k-1)*h)*sin(t*(j-1));
    end
             
    cur = B*u +C0*f0+C1*f1;
    u = Ainv*cur;
    
    U(:, j+1) = u;
              
end


