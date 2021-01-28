format long
N = 32;
M = 32;
V = 1/2;
X = linspace(0, 2*pi, N+1);
X = X(1:N)';
u0 = sin(X);
u_etal = u0;

        
%x = linspace(1, N, N);

%result = solve_diff(N, M, V, u0);




[X_mat, T_mat] = meshgrid(X, linspace(0, 2*pi, M+1));
U_exact = sin(X_mat - V*T_mat);
%norma1 = max(abs(result'-U_exact))

%contour(X_mat, T_mat, result', 'ShowText', 'on')
%contour(X_mat, T_mat, result'-U_exact, 'ShowText', 'on')
n1 = 10;
m1 = 10;

error = zeros(1, 10);
mass_n = zeros(1, 10);




for i = 1:4
    n1 = n1*2;
    m1 = m1*2;
    mass_n(i) = n1;
    X = linspace(0, 2*pi, n1+1);
    X = X(1:n1)';
    
    %u0 = sin(X);
    u0 = ones(n1, 1); % если для Кранка-Николсона, то начальное условие sin(X)

    result = solve_diff(n1, m1, V, u0);
    un = result(:, m1+1);
    
    %[X_mat, T_mat] = meshgrid(X, linspace(0, 2*pi, m1+1));
    %U_exact = sin(X_mat - V*T_mat);
    
    
    error(i) = max(abs(un-sin(X)*cos(2*pi)));
end

loglog(mass_n, error);



xlabel('N')
ylabel('error')
disp ((log(error(4))-log(error(1)))/(log(mass_n(4))-log(mass_n(1))))
    
%title(strcat(strcat('N = ', num2str(N))), strcat('M = ', num2str(M)));


   %if j == 1
   %     cur = B*u +C0*f0+C1*f1+C11*f11;
   %     l = Ainv*cur; % u^1
    
   % else
   %     cur1 = A1*l+B*u+C0*f0+C1*f1+C11*f11;
   %     u = l;
   %     l = Ainv_dop*cur1;
        
   % end
    
    
   % U(:, j+1) = u;