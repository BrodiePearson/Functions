function out=make3dredspec(in,stdev)

% First, generate random white (gaussian) spectrum
out=fft3(randn(size(in)));

% construct frequency vectors and multiply
[nx,ny,nz]=size(in);

dk=1/nx;
dl=1/ny;
dm = 1/nz;

[k,l,m]=meshgrid((-nx/2:nx/2-1)*dl,(-ny/2:ny/2-1)*dk,(-nz/2:nz/2-1)*dm);

filt=1./sqrt(k.^2+l.^2+m.^2);
filt(~isfinite(filt))=0;

out=real(ifftn(out.*fftshift(filt)));

out=stdev.*out./std(out(:));
