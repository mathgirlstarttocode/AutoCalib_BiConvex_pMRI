%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = my_group_dct_basis(xsize,k,p)

% Define W*X and W'*Y for a 2d discrete cosine transform W 
% performed on columnwise for matrix X, where each column 
% of x is vec(x*h_i.')

% Input:
%        xsize = a vector of 2 positive integers
%           k  = a positive integer
% Output:
%        W = struct of 2 fields
%            1) W.times: W*x
%            2) W.trans: W'*y

W.times = @(x) my_group_dct2(x,xsize,k,p); 
W.trans = @(x) my_group_idct2(x,xsize,k,p); 


function Y = my_group_dct2(x,xsize,k,p)

for i = 1:size(x,2)
    xx = reshape(x(:,i),[xsize,k]);
    y=dct(dct(xx,[],1),[],2); %dct is unitary, no need to rescale
    Y(:,i) = y(:);
end

function X = my_group_idct2(y,xsize,k,p)

for i = 1:size(y,2)
    yy = reshape(y(:,i),[xsize,k]);
    x=idct(idct(yy,[],1),[],2);
    X(:,i) = x(:);
end
