%%%%%%%%%%
% Custom l1,2 minimization solver
% solve min 1/2t_k \|x - alpha\|^2 + lambda \|x\|_1,2
% t_k is stepsize at iteration k
% lambda is the regularizer
% alpha = x_{k-1} - t_k*grad(f(x_{k-1}))
% for MRI problem, 
% grad(f(x)) = WA'(AW'x-b)
% W is opts.basis, A is forward operator, b is measurement
%
% Y. Ni, "Auto-Calibration and Biconvex Compressive Sensing with
% Applications to Parallel MRI."
%%%%%%%%%%
function Out = my_group_FISTAf(maxit,A,res,k,b,xs,step,lambda_start,proj,xstart)

p = size(b,2);

if numel(res) == 1
    N = res^2;
else
    N = res(1) * res(2);
end

%
if ~isempty(xstart) %start value
    y = xstart;
else   
    %y = rank1_proj(A.trans(b),res^2*k,p); %warm start
    y = rank1_proj(A.trans(b),N,k*p); %warm start
end

x = y;
t0 = 1; %stepsize k 

func_vl = zeros(maxit,1);
rerr = zeros(maxit,1);
best_x = [];
best_rerr = Inf;

lambda = lambda_start;
for i =1:maxit   

    if mod(i, 10) == 0 || i == 1
        fprintf(1, '[Proposed l1,2 FISTA] iter %4d/%4d | lambda = %.3e\n', ...
                i, maxit, lambda);
        drawnow;
    end
    
    lambda = lambda * 0.9^(i-1); %decreasing regularizer
    %lambda = lambda_start * 0.1^(i-1); %decreasing regularizer
    t = step;
    
    grad = A.trans(A.times(x)-b); %calculate grad 
    
    %check stopping
    %if norm(grad) <1e-10
       % break
    %end
    
    alpha = x - t*grad; %alpha
    alpha_rs = reshape(alpha,N,[]);
    nl = vecnorm(alpha_rs');
    
    scale = (nl- t*lambda) ./ nl;
   
    %%update
    x1 = reshape(x,N,[]);
    
    idx0 = find(nl < t*lambda);
    x1(idx0,:)=0;
    idx1 = find(nl >= t*lambda);
    x1(idx1,:) = scale(idx1).' .* alpha_rs(idx1,:);
    
    %another projection
    if proj == true
        [U,S,V] = svd(x1,'econ'); 
        x1 = U(:,1)*S(1)*V(:,1)';
    end
    
    %check function value
    %func_vl(i) = evaluate_f(t,lambda,x,alpha,res,k);
    func_vl(i) =  evaluate_f2(A,b,lambda,x,res,k,p);
    if ~isempty(xs) %if simulated data, calculate rel_l2 error
        rerr(i) = norm(x - xs,'fro') / norm(xs,'fro');
    end
    %update
    y1 = reshape(x1,N*k,p);
    t1 = (1+sqrt(1+4*t0^2))/2;
    x1 = y1 + (t0 - 1)/t1 * (y1 - y);
    
    %store best so far
    if rerr(i) < best_rerr
        best_x = x1;
        best_rerr = rerr(i);
    end
    
    %
    %if i>1 && rerr(i) > rerr(i-1)
       %break
    %end
    
    %
    x = x1;
    y = y1;
    t0 = t1;
    
end
Out.x = x;
Out.fv = func_vl;
Out.rerr = rerr;
Out.best_x = best_x;
Out.best_rerr = best_rerr;
Out.lambda = lambda_start;
end


function y = evaluate_f2(A,b,lambda,x,res,k,p)
%evaluate function value 
%
y = 1/2 * norm(A.times(x)-b,'fro').^2 + lambda*my_group_12(x,res,k,p);
end

function y = rank1_proj(x,m,n)
x0 = reshape(x,m,n);
[U,S,V] = svd(x0,'econ');
y = U(:,1) * S(1) * V(:,1)'; 
y = reshape(y,size(x));
end