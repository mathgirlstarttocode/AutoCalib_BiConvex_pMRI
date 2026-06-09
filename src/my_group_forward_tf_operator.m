%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = my_group_forward_tf_operator(B,W,picks,perm,res)

% Define A*x and A'*y for a partial DFT matrix A
% Input:
%            n = integer > 0
%        picks = sub-vector of a permutation of 1:n
% Output:
%        A = struct of 2 fields
%            1) A.times: A*x
%            2) A.trans: A'*y

A.times = @(x) pdft2_n2m(scale_by_B(W.trans(x),B),picks,perm,res);
A.trans = @(y)W.times(descale_by_B(pdft2_m2n(y,picks,perm,res),B));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = pdft2_n2m(x,picks,perm,res)

% Calculate y = A*x,
% where A is m x n, and consists of m rows of the 
% n by n discrete-Fourier transform (FFT2) matrix.
% The row indices are stored in picks.

%if size(x,2) == 1
%    n = length(x);
%    x = reshape(x,[sqrt(n),sqrt(n)]);
%end

n = size(x,1);

if nargin < 4
    res_1 = sqrt(n);
    res_2 = sqrt(n);
else
    res_1 = res(1);
    res_2 = res(2);
end

for i = 1:size(x,2)
    xx = x(:,i);
    xx = reshape(xx(perm),[res_1,res_2]);

%[picks_row,picks_cl]= ind2sub([sqrt(n),sqrt(n)],picks);
%tx = fft(x(perm))/sqrt(n);

    tx = fft2(xx)/sqrt(res_1)/sqrt(res_2);
    tx = tx(:);
    y(:,i) = tx(picks);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = pdft2_m2n(y,picks,perm,res)

% Calculate x = A'*y,
% where A is m x n, and consists of m rows of the 
% n by n inverse discrete-Fourier transform (IFFT) 
% matrix. The row indices are stored in picks.

n = length(perm); 
if nargin < 4
    res_1 = sqrt(n);
    res_2 = sqrt(n);
else
    res_1 = res(1);
    res_2 = res(2);
end

for i = 1 : size(y,2)
    tx = zeros(n,1);
    tx(picks) = y(:,i);
    %x = ifft2(reshape(tx,[res,res]))*sqrt(n);
    xx(perm) = ifft2(reshape(tx,[res_1,res_2]))*sqrt(res_1)*sqrt(res_2);
    x(:,i) = xx(:);
end
%x = reshape(x,[res,res]);

function y = scale_by_B(x,B)
[N,k] = size(B);
B_lft = spdiags(B,[0:N:(k-1)*N],N,k*N);
y = B_lft * x ;



function x = descale_by_B(y,B)

[N,k] = size(B);
B_lft = spdiags(conj(B),[0:-N:-(k-1)*N],k*N,N);
x = B_lft * y ;
