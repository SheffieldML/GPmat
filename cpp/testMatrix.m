
% Inverse tests
A = randn(10);
Ainv = inv(A);

save testInv.mat A Ainv

% Cholesky
C = A*A';
U = chol(C);
L = U';

save testCholesky.mat C L U

% GEMM
D = randn(10, 20);
E = randn(20, 10);
F = randn(10);
G = randn(20);
H = randn(10, 20);
alpha = 3.1;
beta = 2.2;

GEMM1 = alpha*D*E + beta*F;
GEMM2 = alpha*D'*E' + beta*G;
GEMM3 = alpha*D*H' + beta*GEMM1;
GEMM4 = alpha*D'*H + beta*GEMM2;

save testGemm.mat D E F G H alpha beta GEMM1 GEMM2 GEMM3 GEMM4

% TRMM
L = randn(15);
L = tril(L);
LU = L-diag(diag(L)) + eye(size(L));
L2 = randn(30);
L2 = tril(L2);
L2U = L2-diag(diag(L2))+eye(size(L2));
U = randn(15);
U = triu(U);
UU = U-diag(diag(U)) + eye(size(U));
U2 = randn(30);
U2 = triu(U2);
U2U = U2-diag(diag(U2))+eye(size(U2));
B = randn(15, 30);
alpha = pi*exp(1);

TRMM1 = L*B*alpha;
TRMM2 = L'*B*alpha;
TRMM3 = B*L2*alpha;
TRMM4 = B*L2'*alpha;
TRMM5 = LU*B*alpha;
TRMM6 = LU'*B*alpha;
TRMM7 = B*L2U*alpha;
TRMM8 = B*L2U'*alpha;
TRMM9 = U*B*alpha;
TRMM10 = U'*B*alpha;
TRMM11 = B*U2*alpha;
TRMM12 = B*U2'*alpha;
TRMM13 = UU*B*alpha;
TRMM14 = UU'*B*alpha;
TRMM15 = B*U2U*alpha;
TRMM16 = B*U2U'*alpha;

save testTrmm.mat L L2 U U2 B alpha TRMM1 TRMM2 ...
    TRMM3 TRMM4 TRMM5 TRMM6 TRMM7 TRMM8 TRMM9 TRMM10 ...
    TRMM11 TRMM12 TRMM13 TRMM14 TRMM15 TRMM16

% TRSM
L = randn(15);
L = tril(L);
LU = L-diag(diag(L)) + eye(size(L));
L2 = randn(30);
L2 = tril(L2);
L2U = L2-diag(diag(L2))+eye(size(L2));
U = randn(15);
U = triu(U);
UU = U-diag(diag(U)) + eye(size(U));
U2 = randn(30);
U2 = triu(U2);
U2U = U2-diag(diag(U2))+eye(size(U2));
B = randn(15, 30);
alpha = pi*exp(1);

TRSM1 = L\(B*alpha);
TRSM2 = L'\(B*alpha);
TRSM3 = (B*alpha)/L2;
TRSM4 = (B*alpha)/L2';
TRSM5 = LU\(B*alpha);
TRSM6 = LU'\(B*alpha);
TRSM7 = (B*alpha)/L2U;
TRSM8 = (B*alpha)/L2U';
TRSM9 = U\(B*alpha);
TRSM10 = U'\(B*alpha);
TRSM11 = (B*alpha)/U2;
TRSM12 = (B*alpha)/U2';
TRSM13 = UU\(B*alpha);
TRSM14 = UU'\(B*alpha);
TRSM15 = (B*alpha)/U2U;
TRSM16 = (B*alpha)/U2U';

save testTrsm.mat L L2 LU L2U U U2 UU U2U B alpha TRSM1 TRSM2 ...
    TRSM3 TRSM4 TRSM5 TRSM6 TRSM7 TRSM8 TRSM9 TRSM10 ...
    TRSM11 TRSM12 TRSM13 TRSM14 TRSM15 TRSM16

% SYRK
A = randn(10, 15);
alpha = pi*exp(1);
beta = pi/2;
C = randn(10);
C = C + C';
D = randn(15);
D = D*D';

SYRK1=alpha*A*A' + beta*C;
SYRK2=alpha*A'*A + beta*D;

save testSyrk.mat A C D alpha beta SYRK1 SYRK2

% AXPY

alpha = 2.13;
i = 7;
j = 10;
k = 6;

A = randn(20, 10);
B = randn(20, 10);

C = randn(10, 20);
AXPY1 = B;
AXPY2 = B;
AXPY3 = B;
AXPY4 = B;
AXPY1(i, :) = A(k, :)*alpha + B(i, :);
AXPY2(i, :) = C(:, j)'*alpha + B(i, :);
AXPY3(:, j) = A(:, k)*alpha + B(:, j);
AXPY4(:, j) = C(i, :)'*alpha + B(:, j);

D = randn(20);
AXPY5 = D + alpha*diag(C(i, :));
AXPY6 = D + alpha*diag(B(:, j));
save testAxpy.mat i j k alpha A B C D AXPY1 AXPY2 AXPY3 ...
    AXPY4 AXPY5 AXPY6

% GEMV
alpha=pi;
beta = exp(1);
i = 4;
j = 8;
k = 6;

A = randn(20, 15);
B = randn(15, 15);
C = randn(15, 20);
D = randn(15, 20);
E = randn(20, 15);;
GEMV1 = B;
GEMV2 = B;
GEMV3 = B;
GEMV4 = B;
GEMV5 = B;
GEMV6 = B;
GEMV7 = B;
GEMV8 = B;

% gemvRowRow
GEMV1(i, :) = alpha*D(k, :)*C' + beta*B(i, :);
GEMV2(i, :) = alpha*D(k, :)*E + beta*B(i, :);
% gemvRowCol
GEMV3(i, :) = alpha*A(:, k)'*C' + beta*B(i, :);
GEMV4(i, :) = alpha*A(:, k)'*E + beta*B(i, :);
% gemvColCol
GEMV5(:, j) = alpha*C*A(:, k) + beta*B(:, j);
GEMV6(:, j) = alpha*E'*A(:, k) + beta*B(:, j);
% gemvColRow
GEMV7(:, j) = alpha*C*D(k, :)' + beta*B(:, j);
GEMV8(:, j) = alpha*E'*D(k, :)' + beta*B(:, j);


save testGemv.mat i j k alpha beta A B C D E GEMV1 GEMV2 GEMV3 ...
    GEMV4 GEMV5 GEMV6 GEMV7 GEMV8

% GER
alpha = pi*pi;
i = 2;
j = 4;
k = 12;

A = randn(15, 17);
B = randn(15, 17);
C = randn(17, 15);
D = randn(17, 15);
x = randn(15, 1);
y = randn(17, 1);

GER1 = A + alpha*x*y';
GER2 = A + alpha*C(i, :)'*B(k, :);
GER3 = A + alpha*C(i, :)'*D(:, j)';
GER4 = A + alpha*B(:, j)*D(:, k)';
GER5 = A + alpha*B(:, j)*B(i, :);


save testGer.mat i j k alpha A B C D x y GER1 GER2 GER3 GER4 GER5

% SYR
alpha = pi*pi;
i = 2;
j = 4;

A = randn(15, 15);
A = A + A';
B = randn(15, 15);
x = randn(15, 1);

SYR1 = A + alpha*x*x';
SYR2 = A + alpha*B(i, :)'*B(i, :);
SYR3 = A + alpha*B(:, j)*B(:, j)';


save testSyr.mat i j alpha A B x SYR1 SYR2 SYR3 

% SYMV
alpha = pi;
beta = exp(1);
i = 4;
j = 8;
k = 6;

A = randn(15, 15);

B = randn(15, 15);
C = randn(15, 15);
C = C + C';
D = randn(15, 15);
SYMV1 = B;
SYMV2 = B;
SYMV3 = B;
SYMV4 = B;
SYMV5 = B;
SYMV6 = B;
SYMV7 = B;
SYMV8 = B;

% symvRowRow
SYMV1(i, :) = alpha*D(k, :)*C + beta*B(i, :);
SYMV2(i, :) = alpha*D(k, :)*C + beta*B(i, :);
% symvRowCol
SYMV3(i, :) = alpha*A(:, k)'*C + beta*B(i, :);
SYMV4(i, :) = alpha*A(:, k)'*C + beta*B(i, :);
% symvColCol
SYMV5(:, j) = alpha*C*A(:, k) + beta*B(:, j);
SYMV6(:, j) = alpha*C*A(:, k) + beta*B(:, j);
% symvColRow
SYMV7(:, j) = alpha*C*D(k, :)' + beta*B(:, j);
SYMV8(:, j) = alpha*C*D(k, :)' + beta*B(:, j);


save testSymv.mat i j k alpha beta A B C D SYMV1 SYMV2 SYMV3 ...
    SYMV4 SYMV5 SYMV6 SYMV7 SYMV8

% SCALE
A = randn(12, 9);
B = randn(8, 11);
C = randn(17, 15);
i = 2;
j = 4;
alpha = pi/exp(1);

SCALE1 = A*alpha;
SCALE2 = B;
SCALE3 = C;
SCALE2(i, :) = B(i, :)*alpha;
SCALE3(:, j) = C(:, j)*alpha;

save testScale.mat i j alpha A B C SCALE1 SCALE2 SCALE3
