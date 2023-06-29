function [s_out,y_out,x_out]=alfft2d(s,xy,kmax,dh_out,varargin)
% ALFFT2D ... 2D anti leakage fast fourier transform for data reconstruction
%
% [s_out,x_out]=alfft2d(s,x,dx_out,kmax,varargin)
%
% ALFFT1D uses Matlab's nufft command to apply non-uniform fft to
% reconstructing the reguly sampled 2D signal from the non-unifomly sampled
% input signal. The spectrum obtained using nufft is a statistical
% representation of the conventional spectrum. We iteratively extract
% the most nk energetic point(s) from the spectrum to subtract their
% contributions from the input signal. Please refer to the fowllowing
% article for more details:
% Xu, Sheng, Yu Zhang, Don Pham, Gilles Lambaré, 2005,
% Antileakage Fourier transform for seismic data regularization:
% GEOPHYSICS, 70, No. 4, V87-V95.
%
% OBlIGATORY INPUT:
% s ... input irregularly sampled signal (single channel)
% xy ... coordinate vectors [ay,ax] for s
% kmax ... targeted maximum wavenumbers [kymax,kxmax]
% dh_out ... desired ouput sampling intervals [dy,dx]
%
% OPTIONAL INPUT:
% Value-Name pairs that determine the reconstruction
% 'tol' ... tolerance (value inbetween 0.05-0.001 is recommended)
% ********** default is 0.01 ***************
% 'maxiter' ... max number of iterations
% ********** default is 200 ***************
% 'nk' ... number of most energetic point(s) in the spectrum to be
% extracted during each iteration. Larger for faster convergence, smaller
% for higher accuracy.
% *********** default n=1 ****************
% 'verbose' ... 1 for activating verbose, 0 for muting verbose
% *********** default =1 *****************
%
% OUTPUT:
% s_out ... output regularly sampled signal.
% y_out ... output regularly sampled coordinate, of size(ny,1);
% x_out ... output regularly sampled coordinate, of size(nx,1);
%
% by Yubing Li, Institute of Acoustics, Chinese Academy of Sciences
% email: liyubing@mail.ioa.ac.cn
% last updated on 21 June 2023
%
% NOTE: By using this software, you are agreeing to the terms detailed in
% this software's Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by Yubing Li, Institute of Acoustics, 
% Chinese Academy of Sciences, Beijing, Chine. The copyright and ownership 
% is held by its 'AUTHOR' (identified above). 
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any
%    purpose with the exception of re-selling or re-distributing the
%    SOFTWARE.
%
% 2) The AUTHOR must be acknowledged in any resulting publications or
%    presentations.
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. The AUTHOR makes no warranties or
%    representation as to its accuracy, completeness, or fitness for
%    any purpose. The AUTHOR are under no obligation to provide support
%    of any kind for this SOFTWARE.
%
% 4) We periodically add, change, improve or update this SOFTWARE without
%    notice. Users may contact the AUTHOR via email to ask for the lastest
%    version.
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

% check input parameters
if(nargin<4)
    disp('input parameters are not complete!')
    return
end
if(~isvec(s))
    disp('expecting a vector for signal s!')
    return
end
if(length(s)*2~=length(xy(:)))
    disp('expecting the length of x to be twice as of s!')
    return
end
if(length(kmax)~=2)
    disp('expecting a kmax length to be 2, i.e. [kymax,kxmax]!')
    return
end
% if(kmax(1)>=0.5/dh_out(1) || kmax(2)>=0.5/dh_out(2))
%     disp('kmax shall be less than half of sampling wavenumber (0.5/dx_out)!')
%     return
% end

% load input parameters
tol=1e-2;
maxiter=200;
nk=1;
verbose=1;
nargs=length(varargin);
for k=1:2:nargs
    arg=varargin{k};
    if(~ischar(arg))
        error(['additional arguments must be (string,value) pairs. Entry'...
            int2str(k) ' expected to be a string'])
    end
    recognized=0;
    if(strcmp(arg,'tol'))
        tol=varargin{k+1};
        recognized=1;
    end
    if(strcmp(arg,'maxiter'))
        maxiter=round(varargin{k+1});
        if(maxiter<1 )
            error('maxiter must be an integer greater than 0');
        end
        recognized=1;
    end
    if(strcmp(arg,'nk'))
        nk=round(varargin{k+1});
        if(nk<1 )
            error('nk must be an integer greater than 0');
        end
        recognized=1;
    end
    if(strcmp(arg,'verbose'))
        verbose=round(varargin{k+1});
        if(verbose~=0 && verbose~=1)
            error('verbose must be 1 or 0');
        end
        recognized=1;
    end
    if(~recognized)
        error([arg ' is not a recognized parameter name'])
    end
end

% scaling coordinate values
xlen=length(xy(:));
yy=xy(1:xlen/2)';
xx=xy(xlen/2+1:end)';
s=s(:);

x_min=min(xx);
x_max=max(xx);

y_min=min(yy);
y_max=max(yy);

nx_out=floor((x_max-x_min)/dh_out(2))+1;
ny_out=floor((y_max-y_min)/dh_out(1))+1;

qfy=0:1/ny_out:(ny_out-1)/ny_out;
qfx=0:1/nx_out:(nx_out-1)/nx_out;

y_in=(yy-y_min)/(y_max-y_min)*(ny_out-1);
y_tmp=0:1:ny_out-1; y_tmp=y_tmp';
y_out=qfy'*(y_max-y_min)+y_min;%+xx;
dy_true=y_out(2)-y_out(1);

x_in=(xx-x_min)/(x_max-x_min)*(nx_out-1);
x_tmp=0:1:nx_out-1; x_tmp=x_tmp';
x_out=qfx'*(x_max-x_min)+x_min;%+xx;
dx_true=x_out(2)-x_out(1);

afy=qfy/max(qfy)/dy_true;
afy=-max(afy)/2:1/ny_out:max(afy)/2;
afx=qfx/max(qfx)/dx_true;
afx=-max(afx)/2:1/nx_out:max(afx)/2;

[ymesh,xmesh]=ndgrid(y_tmp,x_tmp);

iter=1;
s_res=s; s_tmp=zeros(size(s_res));
s_iter=zeros(ny_out,nx_out);
s_out=s_iter;

if(verbose==1)
    figure(1);figuresize(10,10,'inches');
end

if(kmax(1)<=0.5/dh_out(1))
    kxrange=1:round(kmax(2)*dx_true*(nx_out-1));
end
if(kmax(2)<=0.5/dh_out(2))
    kyrange=round(kmax(1)*dy_true*(ny_out-1))+1:ny_out-round(kmax(1)*dy_true*(ny_out-1));
end

% iterative reconstruction
while mean(abs(s_res))>max(abs(s))*tol && iter<=maxiter

    %     for ii=1:length(x_in)
    %         s_tmp(ii)=interp2(x_tmp,y_tmp,s_iter,x_in(ii),y_in(ii));
    %     end
    F=griddedInterpolant(ymesh,xmesh,s_iter);
    s_tmp(:)=F(y_in,x_in);
    s_res=s_res-s_tmp;
    s_out=s_out+s_iter;

    ks=nufftn(s_res,[y_in,x_in],{qfy,qfx});
    ks=reshape(ks,length(qfy),length(qfx));
    ks_amp=abs(ks);
    ks_tmp=zeros(size(ks_amp));
    if(kmax(1)<=0.5/dh_out(1))
        ks_tmp(:,kxrange)=ks_amp(:,kxrange);
    else
        ks_tmp=ks_amp;
    end
    if(kmax(2)<=0.5/dh_out(2))
        ks_tmp(kyrange,:)=0;
    end
    ks_sparse=zeros(size(ks));
    for ipos=1:nk
        [~,pos]=max(ks_tmp(:)); ks_tmp(pos)=0;
        ks_sparse(pos)=ks(pos);
        posy=mod(pos,size(ks,1));
        posx=ceil(pos/ny_out);
        if  posy<size(ks,1)/2
            if posx>1 && posy>1
                ks_sparse(end-posy+2,end-posx+2)=conj(ks(posy,posx));
            end
        else
            if posx>1
                ks_sparse(size(ks,1)-posy+2,end-posx+2)=conj(ks(posy,posx));
            end
        end
    end
    if isreal(s)
        s_iter=real(ifft2(ks_sparse));
    else
        s_iter=ifft2(ks_sparse);
    end

    if(verbose==1)
        figure(1)

        subplot 221
        if isreal(s_out)
            scatter(xx,yy,5,s_res,'filled'); set(gca,'ydir','reverse')
        else
            scatter(xx,yy,5,abs(s_res),'filled'); set(gca,'ydir','reverse')
        end
        xlim([x_min x_max]); ylim([y_min y_max])

        subplot 222
        if isreal(s_out)
            imagesc(x_out,y_out,s_out);
        else
            imagesc(x_out,y_out,abs(s_out));
        end
        xlim([x_min x_max]); ylim([y_min y_max])

        subplot 223
        imagesc(afx,afy,fftshift(ks_amp));
        xlim([max(-kmax(2)*1.5,afx(1)),min(afx(end),kmax(2)*1.5)])
        ylim([max(-kmax(1)*1.5,afy(1)),min(afy(end),kmax(1)*1.5)])

        subplot 224
        imagesc(afx,afy,fftshift(abs(fft2(s_out))))
        xlim([max(-kmax(2)*1.5,afx(1)),min(afx(end),kmax(2)*1.5)])
        ylim([max(-kmax(1)*1.5,afy(1)),min(afy(end),kmax(1)*1.5)])

        disp(['iteration=',num2str(iter)]);

        drawnow
    end

    iter=iter+1;

end

if(verbose==1)
    figure(1);
    subplot 221
    xlabel('X (m)');ylabel('Y (m)')
    title('Residual at sampling points')
    subplot 222
    xlabel('X (m)');ylabel('Y (m)')
    title('Reconstruction')
    subplot 223
    xlabel('Kx (m^{-1})'); ylabel('Ky (m^{-1})')
    title('Residual NUFFT spectrum')
    subplot 224
    xlabel('Kx (m^{-1})'); ylabel('Ky (m^{-1})')
    title('Reconstructed spectrum')
end