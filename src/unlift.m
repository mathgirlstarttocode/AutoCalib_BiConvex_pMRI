function [x_recon,rerr]=unlift(X,im_nl,res,p,k)

%
x_recon = zeros(res,res);

%{
for i = 1: p 
    [U,~,~] = svd(reshape(X(:,i),res^2,k),'econ');
   
    xx = abs(U(:,1));
    yy = dot(xx(:),abs(im_nl(:)))/norm(xx(:))^2;
    x_recon = x_recon+reshape(xx .*yy,res,res);
end
%}

X0 = reshape(X,res^2,k*p);
[U,S,V] = svd(X0,'econ');
x_recon = reshape(U(:,1),size(im_nl));

if norm(x_recon(:)) ~= 0 
    yy = dot(x_recon(:),abs(im_nl(:)))/norm(x_recon(:))^2;
    x_recon = x_recon .* yy;
end

rerr = norm(x_recon - abs(im_nl),'fro')/norm(abs(im_nl),'fro');

end
