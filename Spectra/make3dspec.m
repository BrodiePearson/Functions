function out=make3dspec(in,stdev,n)

out = fftn(randn(size(in)));

[nx,ny,nz] = size(in);

dk = 1/nx;
dl = 1/ny;
dm = 1/nz;

[k,l,m]=meshgrid((-ny/2:ny/2-1)*dl,(-nx/2:nx/2-1)*dk,(-nz/2:nz/2-1)*dm);
%[k,l,m]=meshgrid((-nx/2:nx/2-1)*dk,(-ny/2:ny/2-1)*dl,(-nz/2:nz/2-1)*dm);

filt=sqrt(k.^2 + l.^2 + m.^2).^(n/2);
filt(~isfinite(filt))=0;

out=real(ifftn(out.*fftshift(filt)));

out=stdev.*out./std(out(:));
