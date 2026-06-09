
% This script uses MRIPhantomv0-8 to generate the Shepp-Logan / brain phantom.
% MRIPhantomv0-8 is not included in this repository. Please download it
% separately and set `phantom_root` below to the local path on your machine.

phantom_root = '/path/to/MRIPhantomv0-8';  % user should edit this
addpath(phantom_root);
addpath(fullfile(phantom_root, 'misc'));

% Set project root as the parent folder of this script's folder
this_file = mfilename('fullpath');
[this_dir, ~, ~] = fileparts(this_file);
project_root = fileparts(this_dir);

addpath(project_root);
addpath(fullfile(project_root, 'src'));
addpath(genpath(fullfile(project_root, 'external')));
data_dir = fullfile(project_root, 'data');
output_dir = fullfile(project_root, 'output');


noise_spec = 0.05; %gaussian noise level of 0,0.1%,1%,2.5%,10%
model_violation = true; %allow model violation (due to coil only approximated by Bh)

%set parameters
opts.res = 256;  
opts.tol = 1e-6;
opts.maxit = 600;
opts.print = 1;
opts.noise= noise_spec;
opts.Nb_coil = 8;

coil_opts.res = opts.res;
coil_opts.type = 'biot';
coil_opts.param.FS = 0.28; % FOV width
coil_opts.param.D = 0.17; % Distance center->coil
coil_opts.param.R = 0.09; % radius of the coil
coil_opts.sens.model ='polynomial'; % 'sinusoidal','polynomial';
coil_opts.sens.param = 2; %number of basis to span coil sensitivity
coil_opts.Nb_coils = opts.Nb_coil ;
coil_opts.param.rand = true; %add randomness to the coils
save(fullfile(data_dir, 'phantom', 'opts.mat'), "opts", "coil_opts")


% PHANTOM
DefineBrain;
DefineSL;
leg = {'brain', 'Sheep-Logan'};
im = RasterizePhantom(Brain,opts.res,[1],0);
save(fullfile(data_dir, 'phantom', 'phantom.mat'), 'im');

%generate data
%[h,M,sensitivity]= generate_group_data_v1(opts.res, opts.Nb_coil, coil_opts,im); %by physical law
[im,h,M,sensitivity]= generate_group_data(opts.res,opts.Nb_coil,coil_opts);
save(fullfile(data_dir, 'phantom', 'coil_v2.mat'), "sensitivity", "M", "h")
%generate full k-space
slice_full = fft2(sensitivity .* im)/opts.res;

%generate noise
Signal_Level = rssq(slice_full(:));
Gaussian_noise = randn(size(slice_full));
Noise_Level = rssq(Gaussian_noise(:));
weight = opts.noise*Signal_Level/Noise_Level;
noise = weight * Gaussian_noise;
SNR = Signal_Level/rssq(noise(:));

slice_full = slice_full+ noise;
save(fullfile(data_dir, 'phantom', 'slice_full_v2.mat'), "slice_full", "noise")

%generate full k-space
coil = reshape(M*h, opts.res,opts.res,opts.Nb_coil);
slice_full = fft2(coil .* im)/opts.res;

%generate noise
Signal_Level = rssq(slice_full(:));
Gaussian_noise = randn(size(slice_full));
Noise_Level = rssq(Gaussian_noise(:));
weight = opts.noise*Signal_Level/Noise_Level;
noise = weight * Gaussian_noise;
SNR = Signal_Level/rssq(noise(:));

slice_full = slice_full+ noise;
save(fullfile(data_dir, 'phantom', 'slice_full_v2_noviolation.mat'), "slice_full", "noise")



