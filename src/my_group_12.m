function y = my_group_12(x,res,k,p)
%calculate the mixed 1,2 norm for x
% x = im(:) * h.'; 

if numel(res) == 1
    N = res^2;
else
    N = res(1) * res(2);
end

x = reshape(x,N,k*p);
y = sum(vecnorm(x.'));
end
