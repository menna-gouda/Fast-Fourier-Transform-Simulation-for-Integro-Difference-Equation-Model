function [BICs,Parms,Flags]=fit(y)

alpha=mean(y);d=1;

[parms1, fval1, iflag1]= fminsearch(@(x) nll_bessel(y, x), alpha); 
[parms2, fval2, iflag2]= fminsearch(@(x) nll_typeII(y, x), [alpha d]); 
[parms3, fval3, iflag3]= fminsearch(@(x) nll_typeIII(y, x), [alpha d]); 
[parms4, fval4, iflag4]= fminsearch(@(x) nll_model4(y, x), [alpha d]); 

N=length(y);
BIC1 = 2*fval1 + 1*log(.5*N/pi);
BIC2 = 2*fval2 + 2*log(.5*N/pi);
BIC3 = 2*fval3 + 2*log(.5*N/pi);
BIC4 = 2*fval4 + 2*log(.5*N/pi); 

BICs=[BIC1,BIC2,BIC3,BIC4];
Parms=[parms1,parms2,parms3,parms4];
Flags=[iflag1,iflag2,iflag3,iflag4];

end