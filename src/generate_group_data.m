function [im,h,M,varargout]= generate_group_data(res,p,coil_opts,varargin)
%% PHANTOM
DefineBrain;
DefineSL;
leg = {'brain', 'Sheep-Logan'};

%%
im = RasterizePhantom(Brain,res,[1],0);

%
if ~isfield(coil_opts,'res');        coil_opts.res = res;       end
if ~isfield(coil_opts,'Nb_coils');   coil_opts.Nb_coils = 1;    end  
if ~isfield(coil_opts.param,'FS');   coil_opts.param.FS = 0.28;    end      
if ~isfield(coil_opts.param,'D');    coil_opts.param.D = 0.17;   end      
if ~isfield(coil_opts.param,'R');    coil_opts.param.R = 0.09;   end      
if ~isfield(coil_opts.sens,'model');  coil_opts.sens.model ='polynomial'; end
if ~isfield(coil_opts.sens,'param');  coil_opts.sens.param =2; end
 
%
coil = simulate_sensitivities(coil_opts); %simulate sensitivities


%spanning sensitivity using some basis

NbCoils = size(coil.sensitivity,3);
%support = (im>1e-3);numel(find(support)); % support of the image, which should be assume to be unknown
support = true(size(im));

sensitivity = coil.sensitivity/max(reshape(abs(coil.sensitivity.*repmat(support,[1,1,NbCoils])),1,numel(coil.sensitivity)));
sens.model =coil_opts.sens.model;
sens.param = coil_opts.sens.param; %controls length of h,len(h)=x^2

for i = 1:p 
    if i == 1 
        [sens_sl,nrmse,ser,maxerror,condi,M] = SensFitting(sensitivity(:,:,1),sens.model,sens.param,support);
         h(:,1) = sens_sl.data;
    else 
        [sens_sl,] = SensFitting(sensitivity(:,:,i),sens.model,sens.param,support);
        h(:,i) = sens_sl.data;
    end
end
if nargout > 3
        varargout{1} = sensitivity; % Store sensitivity in the optional output
end
end


