function E = Spectra_from_SF_1D(D,r)

% take derivative of structure function data
dD = diff(D);
dr = diff(r);
dDdr = dD./dr;

n = length(dDdr);

% discrete sine transform 
E = zeros (n,1);
for i = 1:n
    for j = 1:n
      E(i) = E(i) + sin(pi*i*j/n) * dDdr(j);
    end
end

% multiply by coefficients
E = E.*sqrt(1/(2*n));


end