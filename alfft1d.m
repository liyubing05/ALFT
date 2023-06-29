function [s_out,x_out]=alfft1d(s,x,kmax,dx_out,varargin)
% ALFFT1D ... 1D anti leakage fast fourier transform for data reconstruction
%
% [s_out,x_out]=alfft1d(s,x,dx_out,kmax,varargin)
%
% ALFFT1D uses Matlab's nufft command to apply non-uniform fft to 
% reconstructing the reguly sampled 1D signal from the non-unifomly sampled 
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
% x ... coordinate vector for s
% kmax ... targeted maximum wavenumber
% dx_out ... desired ouput sampling interval inbetween [min(x),max(x)]
%
% OPTIONAL INPUT:
% Value-Name pairs that determine the reconstruction
% 'tol' ... tolerance (value inbetween 0.05-0.001 is recommended)
%  ********** default is 0.01 ***************
% 'maxiter' ... max number of iterations
%  ********** default is 200 ***************
% 'nk' ... number of most energetic point(s) in the spectrum to be
% extracted during each iteration. Larger for faster convergence, smaller
% for higher accuracy.
% *********** default n=1 ****************
% 'verbose' ... 1 for activating verbose, 0 for muting verbose
% *********** default =1 *****************
%
% OUTPUT:
% s_out ... output regularly sampled signal.
% x_out ... output regularly sampled coordinate. The same size as s_out.
%
% by Yubing Li, Institute of Acoustics, Chinese Academy of Sciences
% email: liyubing@mail.ioa.ac.cn
% last updated on 14 June 2023
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
if(~isvec(s) || ~isvec(x))
    disp('expecting a vector for signal s or axis x!')
    return
end
if(~isscalar(kmax))
    disp('expecting a scalar for kmax!')
    return
end
if(kmax>=0.5/dx_out)
    disp('kmax shall be less than half of sampling wavenumber (0.5/dx_out)!')
    return
end

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
x=x(:);
s=s(:);

x_min=min(x);
x_max=max(x);

nx_out=floor((x_max-x_min)/dx_out)+1;
qf=0:1/(nx_out-1):1;

x_in=(x-x_min)/(x_max-x_min)*(nx_out-1);
x_tmp=0:1:nx_out-1; x_tmp=x_tmp';
x_out=qf'*(x_max-x_min)+x_min;%+xx; 
dx_true=x_out(2)-x_out(1);

iter=1;
s_res=s;
s_iter=zeros(size(x_tmp));
s_out=s_iter;

if(verbose==1)
    figure(1);figuresize(10,10,'inches');
end

% iterative reconstruction
while mean(abs(s_res))>max(abs(s))*tol && iter<maxiter

    s_tmp=interp1(x_tmp,s_iter,x_in);
    s_res=s_res-s_tmp;
    s_out=s_out+s_iter;

    ks=nufft(s_res,x_in,qf);
    ks_amp=abs(ks);
    ks_tmp=ks_amp(1:round(kmax*dx_true*(nx_out-1)));
    ks_sparse=zeros(size(ks));
    for ipos=1:nk
        [~,pos]=max(ks_tmp); ks_tmp(pos)=0;
        ks_sparse(pos)=ks(pos);
        if pos>1; ks_sparse(end-pos+2)=conj(ks(pos)); end
    end
    s_iter=real(ifft(ks_sparse));

    if(verbose==1)
        figure(1)

        subplot 411
        hold off;scatter(x,s,'*');
        hold on; plot(x_out,s_out,'linewidth',2);
        xlim([x_min x_max])

        subplot 412
        plot(qf/dx_true, abs(ks_sparse),'linewidth',2)
        xlim([0 min(kmax*1.5,max(qf)/dx_true/2)]);

        subplot 413
        hold off;scatter(x,s_res,'*');
        ylim([-max(abs(s))*1.05 max(abs(s))*1.05])
        xlim([x_min x_max])

        subplot 414
        plot(qf/dx_true, abs(ks),'linewidth',2)
        xlim([0 min(kmax*1.5,max(qf)/dx_true/2)]);

        disp(['iteration=',num2str(iter)]);

        drawnow
    end

    iter=iter+1;

end

if(verbose==1)
    figure(1);
    subplot 413
    xlabel('Time (s)');ylabel('Amplitude')
    title('Residual at sampling points')
    subplot 414
    xlabel('Frequency (Hz)'); ylabel('Amplitude')
    title('Non-uniform fourier transform')
    subplot 411
    xlabel('Frequency (Hz)'); ylabel('Amplitude')
    title('Sampled points vs. Reconstruction')
    subplot 412
    xlabel('Time (s)'); ylabel('Amplitude')
    title('Sparse frequency point(s)')
end


