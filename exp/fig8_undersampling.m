%% Phantom robustness experiment: sampling rate vs. noise level
% Fig. 8 of the manuscript:
% Ni, Y., Strohmer, T. Auto-Calibration and Biconvex Compressive Sensing with Applications to Parallel MRI. 
% J Fourier Anal Appl 32, 32 (2026). 
% https://doi.org/10.1007/s00041-025-10223-1
% This script evaluates the proposed lifted l1,2 reconstruction method on a
% simulated phantom under different undersampling rates and Gaussian noise
% levels.

% Set project root as the parent folder of this script's folder
this_file = mfilename('fullpath');
[this_dir, ~, ~] = fileparts(this_file);
project_root = fileparts(this_dir);

addpath(project_root);
addpath(fullfile(project_root, 'src'));
addpath(genpath(fullfile(project_root, 'external')));

data_dir = fullfile(project_root, 'data');
output_dir = fullfile(project_root, 'output','simulation');



%set parameters
opts.res = 256;  % generate 256 by 256 phantom
opts.tol = 1e-4;
opts.maxit = 300;
opts.print = 0;
opts.Nb_coil = 12;


coil_opts.res = opts.res;
coil_opts.type = 'biot';
coil_opts.param.FS = 0.28; % FOV width
coil_opts.param.D = 0.17; % Distance center->coil
coil_opts.param.R = 0.09; % radius of the coil
coil_opts.param.rand = 0;
coil_opts.sens.model ='polynomial'; % 'sinusoidal';
coil_opts.sens.param = 2; %number of basis to span coil sensitivity
coil_opts.Nb_coils = opts.Nb_coil;

%varying reduction factor and noise level
nn = 4; %reduction factor 1/nn*[1:nn]
noise_spec = [0,0.01,0.02,0.05,0.10]; %varying gaussian noise level
mm = length(noise_spec);

%generate data
[im,h,M]= generate_group_data(opts.res,opts.Nb_coil,coil_opts); %by physical law

opts_all =[];
for i = 1:nn 
    for j = 1:mm
opts_all(i,j).opts = opts;
opts_all(i,j).opts.samp_rate = 1/nn * i;
%opts_all(i,j).opts.noise= get_noise_perc(noise_spec(j),opts.res);
opts_all(i,j).opts.noise= noise_spec(j);
    end
end

Rec_Coil = [];
parfor i = 1:nn
    for j = 1:mm
        %solver
        Rec_Coil(i,j).Out_best = lsense(im,h,M,opts_all(i,j).opts,coil_opts,false,[],[],[],[]);
    end
end


plotter(nn,mm,Rec_Coil,opts)

%save(fullfile(output_dir, 'undersampling_simulation_04.mat'))

function plotter(nn,mm,Rec_Coil,opts)
%plotting
f = figure();


ii = 0;

[ha,pos] = tight_subplot(nn, mm, [0.05, 0.05],[0.05 .1]);
titleHandle = sgtitle(sprintf('reconstruction results using %d coils',opts.Nb_coil)); 
%set(titleHandle, 'Position', [0.5, 0.95, 0.5]); % Adjust the position [x, y, width]

for i_sr = 1:nn
    sr = 1/nn * i_sr; % Sampling rate
    for j_nl = 1:mm
        ii = ii + 1;
        axes(ha(ii));
        imshow(abs(Rec_Coil(i_sr,j_nl).Out_best.recon_im.im) , []);
        title(sprintf('rerr %.2f', Rec_Coil(i_sr,j_nl).Out_best.recon_im.rerr));
        % Labeling x-axis and y-axis
        if i_sr == nn
            xlabel(sprintf('NSR: %.2f%%', 100*1/Rec_Coil(i_sr,j_nl).Out_best.SNR));
        end
        if j_nl == 1
            ylabel(sprintf('%.2f%% \n sampling', 100*sr),'Rotation',0);
        end
        
    end
end

end
