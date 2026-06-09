%test on the reconstruction under noise and various sampling masks
%Figure 10 of the paper
%this script is for real brain data, with simulated coil sensitivity maps and noise.
%First generate the k-space data and noise using Simulations/Exp3_generate_brain_data.m, then run this script to perform reconstructions and generate the results.
%After running this code, the reconstruction will be saved to 
%/output/brain/brain_recon_US%d.mat

% The data is adapted from the ENLIVE reproducibility data:08_09_lowrank/head_slice1/data/full.cfl/.hdr
% The original BART-format data are from the ENLIVE repository folder `08_09_lowrank/`, specifically `head_slice1/` and `knee/`.

% H. C. M. Holme, S. Rosenzweig, F. Ong, R. N. Wilke, M. Lustig, and M. Uecker,
% “ENLIVE: An Efficient Nonlinear Method for Calibrationless and Robust Parallel Imaging,”
% Scientific Reports 9:3034, 2019.
% 

clear xs

% Set project root as the parent folder of this script's folder
this_file = mfilename('fullpath');
[this_dir, ~, ~] = fileparts(this_file);
project_root = fileparts(this_dir);

addpath(project_root);
addpath(fullfile(project_root, 'src'));
addpath(genpath(fullfile(project_root, 'external')));

data_dir = fullfile(project_root, 'data');
output_dir = fullfile(project_root, 'output');

%load brain data
data_id = 'brain';
datadir = fullfile(data_dir, 'brain');
outputdir = fullfile(output_dir, 'brain');
data = load(fullfile(datadir,'slice-full.mat'));
data_full = squeeze(data.data); %load full k-space data

%
max_it =300;
num_k = 5;
type = 'polynomial'; %type of basis to expand coil sensitivities
proj = true; 
CTN = false;
if CTN == true
    try
        x_start = xtest;
    end
end

%type = 'low_freq';

%load undersampling mask; undersampling mask can be generated using the NLINV code, from GitHub ().
US = 8.5; %4 times undersampling
sr = 1/US; %sampling rate
if US ~= 1
    pat = load(fullfile(datadir, sprintf('pat_%g.mat', US)));
    %pat_cal = load(strcat(datadir,'pat_', num2str(US), '_cal.mat'));
    pat = squeeze(pat.poisson_mask);
    pat_rearange =ifftshift(pat); %rearange the poisson sampling mask
end

%extract problem dimensions
res_1 = size(data_full,1);
res_2 = size(data_full,2);
p = size(data_full,3);
n = res_1*res_2;

%the BART k-space measurements have shifted freq
for i = 1:32
    shift2 = sqrt(res_1)*sqrt(res_2)*fftshift(ifft2(ifftshift(data_full(:,:,i)))); 
    data_full_rearange(:,:,i) = fft2(shift2)/(sqrt(res_1)*sqrt(res_2)); %original k-space has zero-freq at center
end

support = (abs(shift2)>1e-3);numel(find(support));


%undersampling using a mask
if sr == 1
    ksf = data_full_rearange;%fully sampled k-space data
else
    ksf = data_full_rearange .* pat_rearange; %measurements
end

b = reshape(ksf, [], p);%reshape the ksf columnwise? 

%subsampling, need the masking
perm = 1:n;
if sr == 1
    picks = 1:n;
else
    mask_indexes = pat_rearange(:); %reshape
    picks = find(mask_indexes ~=0); %get the index where is nonzero
end

b = double(b(picks,:));

%parameters and options
res_1 =size(ksf,1);
res_2 =size(ksf,2);
res = res_1;
p = size(ksf,3);
xsize = [res_1 res_2];
n = res_1*res_2;


%%%%%%%%define forward operator

%define forward operator
Mtype = {'my_group_forward_tf'};

%B
if strcmp(type, 'polynomial')
    [h,M]= generate_coil_basis(xsize,p,num_k); %only need the polynomial basis M and size
    k = size(h,1);
    M_nl = vecnorm(M); %column normalization
    M_new = M./M_nl;
    B = M_new;
elseif strcmp(type, 'low_freq')
    M = generate_FFT2(res,num_k);
    k = size(M,2) ; 
    M_nl = vecnorm(M); %column normalization
    M_new = M./M_nl;
    B = M_new;

else
    k = num_k;
    a = 240;
    b = 40;
    W_mat = generate_W(res,k,a,b);
    W_nl = vecnorm(W_mat); %column normalization
    B = W_mat./W_nl;
  
end

%transform sparsity W
W = my_group_dct_basis(xsize,k,p);

%A
op_A = eval(['@' Mtype{1} '_operator']);
A = feval(op_A,B,W,picks,perm,xsize);


%calculate the optimal stepsize for ISTA/FISTA
N = size(M,1);
B_lft = spdiags(M_new,[0:N:(k-1)*N],N,k*N); 
L = svds(B_lft,3);
step = 1/(L(1)^2);
%step =2.5208e+03;
%step = 1;

%initialization
if CTN == false 
    for i = 1:p
    coil_img(:,:,i) = sqrt(res_1)*sqrt(res_2)*fftshift(ifft2(ksf(:,:,i))); 
    end
    img_avg = rssq(coil_img,3); %average coil image
    xs0 = repmat(img_avg(:),k,p); %estimate lifted matrix with uniform coil
    xs0 =  W.times(xs0);
else
    xs0 = x_start;
end

%estimate lambda
rd = 1/2*norm(A.times(A.trans(b))-b,'fro')^2;
lambda = 1/2 * rd /(my_group_12(A.trans(b),xsize,k,p)+my_group_12(xs0,xsize,k,p)); 
%lambda =  6.3930e-06;%US4
%lambda = 1.0086e-07; %US7

results = struct('lambda', {}, 'f', {});
% Estimate regularization parameter from the data residual and initial point.
for iii = -3:-3
    %lambda0 = lambda *10^(1-iii);
    lambda0 = lambda;
    fprintf('lambda %i',lambda0);
    Out.f = my_group_FISTAf(max_it,A,xsize,k,b,[],step,lambda0,proj,double(xs0));
    results(iii + 4).lambda = lambda0; % Index adjusted to be non-negative
    results(iii + 4).f = Out.f;
end


for iii = -3:-3
    %norm(W.trans(x) - xs) / norm(xs);
    xtest = results(iii + 4).f.x;
    residue = norm(A.times(xtest)-b,'fro'); %

    %xtest = Out(2).f.x;
    X = W.trans(xtest);%transform back
    X0 = W.trans(A.trans(b));

    x_recon_percoil = zeros(res_1,res_1,p);
    for pp = 1:p
        [U,S,V] = svd(reshape(X(:,pp),res_1^2,k),'econ');
        x_recon_percoil(:,:,pp) = reshape(abs(U(:,1)),res_1,res_1);
    end
    x_recon_all(:,:,iii+4) = rssq(x_recon_percoil,3);
end

%plotting
% Select reconstruction index.
% In the original code, iii = -3 gives recon_idx = 1.
recon_idx = 1;

recon = abs(x_recon_all(:,:,recon_idx));
lambda_val = results(recon_idx).lambda;

% Sampling information
sampling_rate = nnz(pat) / numel(pat);

figure('Color','w', 'Position', [200, 200, 1100, 450]);

tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

% Reconstruction
nexttile;
imshow(recon, []);
axis image off;
colormap gray;
title(sprintf('Proposed reconstruction\nUS = %g, k = %d, \\lambda = %.2e', ...
      US, k, lambda_val), ...
      'FontWeight','normal');

% Sampling mask
nexttile;
imshow(abs(pat), []);
axis image off;
colormap gray;
title(sprintf('Sampling mask\n%.2f%% sampled, %d iterations', ...
      100*sampling_rate, max_it), ...
      'FontWeight','normal');

sgtitle('Lifted l_{1,2} reconstruction and undersampling pattern', ...
        'FontWeight','bold');



% ------------------------------------------------------------
% Save selected reconstruction
% ------------------------------------------------------------
% 
% recon = abs(x_recon_all(:,:,recon_idx));
% 
% if ~exist(outputdir, 'dir')
%     mkdir(outputdir);
% end
% 
% filename = fullfile(outputdir, sprintf('my_recon_US%g.mat', US));
% save(filename, 'recon');
% 
% fprintf('Saved reconstruction to: %s\n', filename);