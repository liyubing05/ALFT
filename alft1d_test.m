% Demo for alft1d
% Yubing Li, Institute of Acoustics, Chinese Academy of Sciences
% email: liyubing@mail.ioa.ac.cn
% last updated on 14 June 2023

%%
clear; close all;

npts=120;
maxval=1600;
fmax=1.2*npts/maxval;
dt=1;

tt=0:1600;
ff=0:1/max(tt):1/(tt(2)-tt(1)); ff=ff';

% case 1 - summation of several harmonic signals
npos=1;
aa=tukeywin(maxval+1,0.7);
Sc = 2*sin(0.012*pi*tt)'.*aa(tt+1);
for ia=1:5
    aa=tukeywin(maxval+1,rand/2+0.2);
    Sc=Sc+2*rand*sin(3*rand*0.03*pi*tt)'.*aa(tt+1);
end
SS = awgn(Sc,10);

% % case 2 - convolution btw random pulses and a band limited wavelet
% npos=1;
% Sc=zeros(size(tt))';
% for ia=1:35
%     Sc(randi(1400,1)+100)=rand([1,1])*40-20;
% end
% Sc=butterfilter(Sc,tt'/1600,'fmin',3,'fmax',60,'phase',0);
% SS=awgn(Sc,20);

%%
[S,tidx]=datasample(SS,npts,'Replace',false);
t=tt(tidx);

figure(4);figuresize(12,3,'inches')
hold off;
plot(tt,SS,'k-.','linewidth',1);
hold on;
plot(tt,Sc,'r','linewidth',1);
scatter(t,S,'b*');
xlabel('Time (s)'); ylabel('Amplitude')
legend('Noisy signal','Clean sigal','Samplings')

%%
tic;
[s_out,t_out]=alfft1d(S,t,fmax,dt,'nk',npos,'maxiter',2000,'tol',1e-3,'verbose',0);
toc;

%%
% stop
figure(2);figuresize(12,5,'inches')
subplot(211)
hold off;
plot(tt,SS,'k','linewidth',1);
hold on;
plot(tt,Sc,'b','linewidth',2);
plot(t_out,s_out,'r','linewidth',2)
xlabel('Time (s)'); ylabel('Amplitude')
legend('Noisy signal','Clean sigal','Reconstruction')

subplot(212)
hold off;
plot(ff,abs(fft(SS)),'k','linewidth',1)
tt=0:1600;
hold on;
plot(ff,abs(fft(Sc)),'b','linewidth',2)
plot(ff,abs(nufft(s_out,t_out,tt/(max(tt)+tt(2)-tt(1)))),'r','linewidth',2)
legend('Noisy signal','Clean sigal','Reconstruction')
xlim([0 fmax])
xlabel('Frequency (Hz)'); ylabel('Amplitude')
