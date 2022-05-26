function [c_n]=fun_FFT(fun,mm,x)

%% Prepare inputs
if size(x,1) > size(x,2); x = x.'; end
if size(mm,1) > size(mm,2); mm = mm.'; end

%% OLD BAD VERSION
% % [z,w] = MacGaussQuad1D(10);
% % xx_i=.5*(x(end)-x(1))*z+(x(end)+x(1))*.5;
% % c_n=zeros(length(mm),1);
% % for ii=1:length(mm)
% %     f_int=fun.*exp(-1i*mm(ii)*x);
% %     ff_i=interp1(x,f_int,xx_i);
% %     Int=.5*(x(end)-x(1))*w*ff_i';
% %     c_n(ii)=Int/(2*pi);
% % end


%% NEW GOOD VERSION
% % nGauss=10;
% % [z,w] = MacGaussQuad1D(nGauss);
% % x=interp1(1:length(x),x,linspace(1,length(x),100));
% % RESOLUTION=20;
% % xx=interp1(1:length(x),x,linspace(1,length(x),RESOLUTION*numel(x)+1));
% % fun_xx=interp1(x,fun,xx);
% % exp_xx=interp1(x,exp(-1i*mm'*x)',xx)';
% % f_int=repmat(fun_xx,length(mm),1).*exp_xx;
% % 
% % c_n=zeros(length(mm),1);
% % nP=nGauss*RESOLUTION;
% % for ii=1:length(mm)
% %     f_int_jj=f_int(ii,:);
% %     Int=0;
% %     for jj=1:nGauss
% %         ind=(jj-1)*nP+1:jj*nP;
% %         xx_jj=xx(ind);
% %         xx_jj_gauss=.5*(xx_jj(end)-xx_jj(1))*z+(xx_jj(end)+xx_jj(1))*.5;
% %         ff_i=interp1(xx,f_int_jj,xx_jj_gauss);
% %         Int=Int+.5*(xx_jj(end)-xx_jj(1))*w*ff_i';
% %     end
% %     c_n(ii)=Int/(2*pi);
% % end

%% BEST VERSION
[z,w] = MacGaussQuad1D(10);
x0 = (x(1:end-1)+x(2:end))*0.5;
h2 = diff(x)*0.5;
xx = z'*h2 + ones(size(z'))*x0;
wh = w'*h2;
xx = xx(:);  wh=wh(:);

% % c_n=zeros(length(mm),1);
% % for ii=1:length(mm)
% %     expmt = exp(-1i*xx*mm(ii))/(2*pi);
% %     yy    = interp1(x,fun,xx').*wh';
% %     c_n(ii)  = yy*expmt;
% % end

yy    = interp1(x,fun,xx').*wh';
expmt = exp(-1i*xx*mm)/(2*pi);
c_n  = (yy*expmt).';


end



function [x,w] = MacGaussQuad1D(n)

%% Compute rule for 1D for interval [-1,1], with the order = n

eps = 2.0e-15;
xstart = -1.0;  xstop = 1.0;

x = zeros(1,n);
w = zeros(1,n);

m  = round((n+1)/2);
xm = 0.5*(xstop+xstart);
xl = 0.5*(xstop-xstart);

if n>1
    for j1=1:m
        z = cos(pi*(j1-0.25)/(n+0.5));
        z1 = 2*z;
        while (abs(z-z1)>eps)
            p1=1.0;
            p2=0.0;
            for j2=1:n
                p3=p2;
                p2=p1;
                p1=((2.0*j2-1.0)*z*p2-(j2-1.0)*p3)/j2;
            end
            pp=n*(z*p1-p2)/(z*z-1);
            z1=z;
            z=z1-p1/pp;
        end
        k=n+1-j1;
        x(j1)=xm-xl*z;
        x(k)=xm+xl*z;
        w(j1)=2*xl/((1-z*z)*pp*pp);
        w(k)=w(j1);
    end
end

if n==1
    x(1) = 0.0;
    w(1) = 2.0;
end

%disp(['   GaussQuad1D: x=[' num2str(x) '], w=[',num2str(w) ']'])
end

