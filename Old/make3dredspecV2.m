function out=make3dredspecV2(in,stdev)

% First, generate random white (gaussian) spectrum
out = fftn(randn(size(in)));

% construct frequency vectors and multiply
[nx,ny,nz] = size(in);

dk = 1/(nx-1);
dl = 1/(ny-1);
dm = 1/(nz-1);

nxhalf = (nx-1)/2;
nyhalf = (ny-1)/2;
nzhalf = (nz-1)/2;

[k,l,m]=meshgrid((-nxhalf:nxhalf)*dk,(-nyhalf:nyhalf)*dl,(-nzhalf:nzhalf)*dm);

filt = 1./sqrt(k.^2+l.^2+m.^2);
filt(~isfinite(filt))=0;

out=real(ifftn(out.*fftshift(filt)));

out=stdev.*out./std(out(:));