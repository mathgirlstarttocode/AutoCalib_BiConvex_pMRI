function [res,tol,maxit,print,noise,samp_rate,Nb_coil] = get_options(opts)

%get parameters for algorithm
res = 256;
tol = 1e-4;
maxit = 9999;
print = 0;
noise = 0; %noise level
samp_rate = 1; %sampling rate
Nb_coil = 1;

if isfield(opts,'res');   res = opts.res;    end  
if isfield(opts,'tol');   tol = opts.tol;    end            
if isfield(opts,'maxit'); maxit = opts.maxit; end
if isfield(opts,'print'); print = opts.print; end        
if isfield(opts,'noise'); noise = opts.noise; end
if isfield(opts,'samp_rate'); samp_rate = opts.samp_rate; end        
if isfield(opts,'Nb_coil'); Nb_coil = opts.Nb_coil; end

end
