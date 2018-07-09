 function sqerror = calc_sqerror(x0,FvFm_exp)
y0 = x0(1:35);
k = x0(36:end);
FvFm_sim = Fluorescence(k,y0);
errors = FvFm_exp - FvFm_sim;
errors = reshape(errors,[],1);
sqerror = errors'*errors;

 end