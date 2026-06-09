%%%%%%%%%%
% Generate simulated, undersampled k-space measurement.
% Then solve the reconstruction using the lifted l12 minimization 
% Customer lifted $l_{1,2}$ solver using my_group_forward_tf_operator
%%%%%%%%%%


function Out_best = lsense(im,h,M,opts,coil_opts,model_violation,sensitivity,pat,noise,warmstart)
%function Out_best = lsense(opts,coil_opts)

if isempty(sensitivity)  % Check if the sixth argument is provided
        model_violation = false; % Default value
end

%get parameters and options
[res,tol,maxit,print,noise_level,samp_rate,num_coil] = get_options(opts);

%generate data
%[im,h,M]= generate_group_data(res,num_coil,coil_opts); %by physical law
im_nl = im;

%
xsize = [res res];
n = prod(xsize);
k = size(h,1);


%plot
if model_violation == false
    sens = M*h;
else
    sens = reshape(sensitivity,res^2,[]);
end

p = num_coil;

if print == 1
    figure
    sgtitle('complex coil sensitivities')
    for i=1:p
        subplot(4,p,i)
        dcolor(reshape(sens(:,i),res,res)) 
        subplot(4,p,i+p)
        dcolor(im .* reshape(sens(:,i),res,res)) 
        subplot(4,p,i+2*p)
        imshow(abs(reshape(sens(:,i),res,res)),[])
        subplot(4,p,i+3*p)
        imshow(im .* abs(reshape(sens(:,i),res,res)),[])
    end
end

avg_sens = reshape(sum(sens,2),res,res);



%normalization
M_nl = vecnorm(M); %basis has normalized column
B = M./M_nl;
h_nl = h.*M_nl';


%transform sparsity
W = my_group_dct_basis(xsize,k,num_coil);

%Lifted Data
for i = 1:num_coil
    XS = im_nl(:) * h_nl(:,i).';
    xs(:,i) = XS(:);
end
xs = W.times(xs); 

%%%%%%%%%%%%random sample
m = max(floor(samp_rate*n),1);

if isempty(pat) == 0
    pp = fftshift(squeeze(pat));
    picks = find(pp == 1);
    m = length(picks);
else
    pp = randperm(n);
    picks = sort(pp(1:m)); 
    bw = 5; %pick center
    pp = fftshift(reshape(1:n,res,res));
    pick0 = pp(floor(res/2)-bw:floor(res/2)+bw,floor(res/2)-bw:floor(res/2)+bw);
    picks(1:numel(pick0)) = pick0(:);
    picks = sort(picks);
end

%perm = randperm(n);
perm = 1:n;
mask = zeros(res,res);
mask(picks) = 1;

%%%%%%%%define forward operator
Mtype = {'my_group_forward_tf'};

op_A = eval(['@' Mtype{1} '_operator']);
A = feval(op_A,B,W,picks,perm,[res,res]);


%generate measurement
%same logic as data_generator. For real data, replace the block with
%measurements
%b0 = A.times(xs) + noise; 
if model_violation == false
    b0 = A.times(xs);
else
    bb0 = sensitivity .* im_nl; 
    for i = 1:p
        tx = reshape(fft2(bb0(:, :, i))/res,[],1); % FFT for each slice
        b0(:,i) = tx(picks);
    end
end

if isempty(noise)
    %generate noise
    Signal_Level = rssq(b0(:));
    Gaussian_noise = randn(m,p);
    Noise_Level = rssq(Gaussian_noise(:));
    weight = noise_level*Signal_Level/Noise_Level;
    noise = weight * Gaussian_noise;
    SNR = rssq(b0(:))/rssq(noise(:));
    b = b0 + noise;
else
    noise = reshape(noise,[],p);
    noise = noise(picks,:);
    SNR = rssq(b0(:))/rssq(noise(:));
    b = b0 + noise;
end  

if isempty(warmstart)
    xstart = [];
else
    coil_img = zeros(res,res,p);
    for i = 1:p
        bb = zeros(res^2,1);
        bb(picks) = b(:,i);
        coil_img(:,:,i) = sqrt(res)*sqrt(res)*ifft2(reshape(bb,res,res)); % Zero-filled inverse FFT warm start for each coil.
    end
    img_avg = rssq(coil_img,3); %average coil image
    % figure(2); clf; %optionally plot the averaged coil image
    % imshow(abs(img_avg), []);
    % axis image off;
    % colormap gray;
    % colorbar;
    % title(sprintf('Warm start: RSS zero-filled IFFT, %d coils, %.2f%% sampled', ...
    %       p, 100*numel(picks)/(res^2)));
    y = repmat(img_avg(:),1,p); %estimate lifted matrix with uniform coil, h= [1,0,...,0]
    xstart =  W.times(padarray(y,res*res*(k-1),0,'post'));
end

%calculate the optimal stepsize for ISTA/FISTA
step = get_step(B,k);
%step =2.5208e-03;

%estimate lambda (parameter for regularizer) for FISTA
lambda = get_lambda(A,b,xs,res,k,num_coil);
%lambda = 4.77314805809719e-05; 
%lambda = 1e-07;
%lambda = 1e-07;

%iterations
j=0;
if noise_level == 0 
    for i = 1:1 %solving FISTA with different lambda
       j = j+1;
       if print == 1
           fprintf('Solving proposed lifted l1,2 reconstruction with lambda %.3e\n',lambda*10^(i));
       end
        Out(j) = my_group_FISTAf(maxit,A,res,k,b,xs,step,lambda*10^(i),true,xstart);
    end
else
    %for i = -4:1 %try more lambda when there is noise
    for i = 0
        j = j+1;
        if print == 1
            fprintf('Solving proposed lifted l1,2 reconstruction with lambda %.3e\n',lambda*10^(i));
        end
        Out(j) = my_group_FISTAf(maxit,A,res,k,b,xs,step,lambda*10^(i),true,xstart);
    end   
end

%get the best recon
Out_best = get_best(Out);
X = Out_best.best_x;
X = W.trans(X);%transform back
X0 = W.trans(A.trans(b));%starting point


%undo the lifting to get back the image
[x_recon,rerr_recon]=unlift(X,im_nl,res,p,k);
[x0,rerr_start]=unlift(X0,im_nl,res,p,k);


Out_best.recon_im.im = x_recon; %store the reconstructed image
Out_best.recon_im.rerr = rerr_recon;
Out_best.org_im = im_nl;
Out_best.SNR = SNR;
Out_best.data = b;
Out_best.samp_pattern = mask;
Out_best.start = x0;
end


function Out_best = get_best(Out)

for i = 1:size(Out,2)
    best_rerr(i) = Out(i).best_rerr;
end

idx = min(find(best_rerr == min(best_rerr)));
Out_best = Out(idx);

end

function step = get_step(B,k)
    
N = size(B,1);
B_lft = spdiags(B,[0:N:(k-1)*N],N,k*N); 
%L = svds(B_lft);
L = svds(B_lft,3,'largest','MaxIterations',500);
step = 1/(L(1)^2);

end

function lambda = get_lambda(A,b,xs,res,k,p)
    
rd = 1/2*norm(A.times(A.trans(b))-b,'fro')^2;
lambda = 1/2 * rd /(my_group_12(A.trans(b),res,k,p)+my_group_12(xs,res,k,p)); 
end


