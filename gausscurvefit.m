function gc = gausscurvefit(x,mu,sigma)

gc=1/sqrt(2*pi)/sigma*exp(-(x-mu).^2/2/sigma/sigma);

end