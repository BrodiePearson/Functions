% function out=make2dredspec(in,stdev)
clear all; clc;close all;
in = NaN([100,100]); stdev = 1;
% First, generate random white (gaussian) spectrum
spec=fft2(randn(size(in)));

% construct frequency vectors and multiply
[nx,ny]=size(in);

dk=1/nx;
dl=1/ny;

[k,l]=meshgrid((-ny/2:ny/2-1)*dl,(-nx/2:nx/2-1)*dk);

filt=1./sqrt(k.^2+l.^2);
filt(~isfinite(filt))=0;

outreal=real(ifft2(spec.*fftshift(filt)));

outreal= (outreal-mean(outreal(:)))./std(outreal(:));

[psd kxx] = pwelch(outreal);
psdxx = nanmean(psd,2);

loglog(kxx,psdxx)
hold on
loglog(kxx(10:end),kxx(10:end).^-2)

%move mean
movmean(psd,.01)

title('Spectra -2')
xlabel('k')

% structure function values
SF= NaN(size(outreal));

for cc = 1:size(SF,2)
    for ii = 1:size(SF,1)
        SF(ii,cc) = nanmean((outreal(:,cc)-circshift(outreal(:,cc),ii)).^2);
    end
end

SFav = nanmean(SF,2);

figure
loglog(SFav)
hold on
% loglog(1:15,0.5*[1:15])

title('Structure Function')
xlabel('r')


