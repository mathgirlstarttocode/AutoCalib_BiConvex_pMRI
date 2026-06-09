function [h,M]= generate_coil_basis(res,p,num_k)
%%
coil.Nb_coils = p;
coil.res = res;
coil.type = 'biot';
coil.param.FS = 0.28; % FOV width
coil.param.D = 0.17; % Distance center->coil
coil.param.R = 0.09; % radius of the coil
coil.param.rand = 0;
coil = simulate_sensitivities(coil);

NbCoils = size(coil.sensitivity,3);
%support = (im>1e-3);numel(find(support)); % support of the image, which should be assume to be unknown
support = true(res);
sensitivity = coil.sensitivity/max(reshape(abs(coil.sensitivity.*repmat(support,[1,1,NbCoils])),1,numel(coil.sensitivity)));

%sens.model = 'sinusoidal';
sens.model ='polynomial';
sens.param = num_k; %controls length of h,len(h)=x^2

for i = 1:p 
    if i == 1 
        [sens_sl,nrmse,ser,maxerror,condi,M] = SensFitting(sensitivity(:,:,1),sens.model,sens.param,support);
         h(:,1) = sens_sl.data;
    else 
        [sens_sl,] = SensFitting(sensitivity(:,:,i),sens.model,sens.param,support);
        h(:,i) = sens_sl.data;
    end
end
end
