% Demo for alft2d
% Yubing Li, Institute of Acoustics, Chinese Academy of Sciences
% email: liyubing@mail.ioa.ac.cn
% last updated on 22 June 2023

%%
clear; close all;

npts=200;
maxval_x=200;
maxval_y=210;
% kxmax=0.25*sqrt(npts)/maxval_x;
% kymax=0.25*sqrt(npts)/maxval_y;
kxmax=0.025;
kymax=0.025;
dx=1;
dy=1;

xx=0:maxval_x; ax=xx/maxval_x;
yy=0:maxval_y; ay=yy/maxval_y;
[xmesh,ymesh]=meshgrid(xx,yy);

% % case 1 - peaks
npos=4;
Sc=peaks(201); Sc(201:maxval_y+1,:)=0;
SS=awgn(Sc,5);

%%
[S,idx]=datasample(SS(:),npts,'Replace',false);
% load dd.mat;
x=xmesh(idx)';
y=ymesh(idx)';

figure(4);figuresize(12,3,'inches')
subplot(131)
imagesc(xx,yy,Sc);
title('Clean image')
xlabel('X (m)'); ylabel('Y (m)')
subplot(132)
imagesc(xx,yy,SS);
title('Noisy image')
xlabel('X (m)'); ylabel('Y (m)')
subplot(133)
scatter(x,y,5,S,'filled'); set(gca,'ydir','reverse')
xlim([xx(1) xx(end)])
ylim([yy(1) yy(end)])
xlabel('X (m)'); ylabel('Y (m)')
title('Sampled image')

%%
tic;
[s_out,y_out,x_out]=alfft2d(S,[y,x],[kymax,kxmax],[dy,dx],'nk',npos,'maxiter',5000,'tol',1e-3,'verbose',0);
toc;
stop
%%
ny_out=length(y_out);
nx_out=length(x_out);
qfy=0:1/ny_out:(ny_out-1)/ny_out;
afy=-max(qfy)/2:1/ny_out:max(qfy)/2;
qfx=0:1/nx_out:(nx_out-1)/nx_out;
afx=-max(qfx)/2:1/nx_out:max(qfx)/2;

figure(2);figuresize(12,10,'inches')
subplot(221)
imagesc(xx,yy,SS);
title('Noisy image')
xlabel('X (m)'); ylabel('Y (m)')

subplot(222)
imagesc(x_out,y_out,s_out);
title('Reconstructed image')
xlabel('X (m)'); ylabel('Y (m)')

subplot(223)
imagesc(afx,afy,fftshift(abs(fft2(SS))));
title('Spectrum for noisy image')
xlim([-kymax*1.5,kymax*1.5])
ylim([-kxmax*1.5,kxmax*1.5])
xlabel('Kx (m^{-1})'); ylabel('Ky (m^{-1})');

subplot(224)
imagesc(afx,afy,fftshift(abs(fft2(s_out))))
title('Spectrum for reconstructed image')
xlim([-kymax*1.5,kymax*1.5])
ylim([-kxmax*1.5,kxmax*1.5])
xlabel('Kx (m^{-1})'); ylabel('Ky (m^{-1})');

