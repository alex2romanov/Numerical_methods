N1 = [5, 40];
M1 = [5, 40];
D = 1;
Z = zeros(length(N1), length(M1));
T = zeros(length(N1), length(M1));

for nn = 1:length(N1)
    for mm = 1:length(M1)
        
        tic % timer
        N = N1(nn);
        M = M1(mm);
        
            
        a0 = 2*(6*D*t/(h*h)-1);
        a1 = 2*(6*D*t/(h*h)+1);
        b0 = -4*(6*D*t/(h*h)+5);
        b1 = -4*(6*D*t/(h*h)-5);
        p0 = -t;
        p1 = -t;
        q0 = -10*t;
        q1 = -10*t;
          

        
%         a0 = 1/6-V*t/4/h;
%         b0 = 2/3;
%         c0 = 1/6+V*t/4/h;
%         a1 = -1/6-V*t/4/h;
%         b1 = -2/3;
%         c1 = -1/6+V*t/4/h;
        

        
        
        %u0 = zeros(N, 1);
        %u_etal = zeros(N, 1);
        
        X = linspace(0, 2*pi, N+1);
        X = X(1:N)';
        u0 = sin(X);
        u_etal = u0;

        
        x = linspace(1, N, N);
        
        u = u0;
        result = zeros(N, M+1);
        result(:, 1) = u0;
        
        
        %for j = 1:M
            
        %    cur = B*u;
        %    u = A\cur;
        %    result(:, j+1) = u;
  
        %end
        
        [X_mat, T_mat] = meshgrid(X, linspace(0, 2*pi, M+1));
        contour(X_mat, T_mat, result')
        
        %if ((N==20) && (M == 40))
         %   hold on
         %   plot(x, abs(u-u_etal))

            
        %end
        
        norma1 = max(abs(u0-u_etal));
        Z(nn, mm) = norma1;
        T(nn, mm) = toc;
        
    end
end
%contour(N1,  M1, T);


