function err=nll_model4(r,x)

alpha=x(1);
d=x(2);

[p,c1]=model4(x);

m=alpha./r;
A=(m.^(2/3)).*nthroot(2+2*sqrt(1+((8/27).*m.^2)),3);
B=(m.^(2/3)).*nthroot(2-2*sqrt(1+((8/27).*m.^2)),3);
t=1./(A+B);

err=length(r)*(2*log(alpha)-log(c1))+sum(-log(t)+(d+2).*log(t+1)+0.5.*log(1+((r./alpha).^2)./(2.*(t.^3)))+(d+1).*(0.5.*(t.^2)-t+(1./(4.*t)).*((r/alpha).^2)));

end