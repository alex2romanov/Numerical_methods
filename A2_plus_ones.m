syms nu;
%nu = 1;
N = 4;
A2 = diag(ones(N, 1)*-4*nu)+diag(ones(N-1,1)*(2+2*nu),1)+diag(ones(N-1,1)*(-2+2*nu),-1);
A2(1, N) = -2+2*nu; A2(N, 1) = 2+2*nu;
A2(1,:)= ones(N, 1);
%fplot(det(A2), [-1 1])