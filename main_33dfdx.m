D = 1;
n1 = 24;
m1 = 24;
X = linspace(0, 2*pi, n1+1);
X = X(1:n1)';
u0 = sin(X)*cos(0);

%u0 = ones(n1, 1);

nu = 2;
U = diff_krank(n1, m1, D, u0, nu);