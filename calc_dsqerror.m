function dsqerror = calc_dsqerror(x0,FvFm_exp,errors)

y0 = x0(1:35);
k = x0(36:end);
JAC = zeros(length(errors),length(x0));
for ivar = 1:length(x0)
    h = zeros(length(x0),1);
    h(ivar) = sqrt(eps);
    xtemp = x0+h;
    ktemp = xtemp(1:35);
    y0temp = xtemp(36:end);
    FvFm_temp = Fluorescence(ktemp,y0temp);
    errors_temp = FvFm_exp - FvFm_temp;
    JAC(ivar,:) = reshape((errors_temp-errors)./sqrt(eps),[],1);
end
dsqerror = 2*JAC*errors;


    
    
    
    
