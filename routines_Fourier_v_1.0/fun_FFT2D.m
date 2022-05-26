function [c_n]=fun_FFT2D(fun_xy,nn,mm,x,y)
% 2D FOURIER SERIES

[X,Y]=meshgrid(x,y);

[z,w] = MacGaussQuad1D(10);
x0 = (x(1:end-1)+x(2:end))*0.5;
hx2 = diff(x)*0.5;
xx = z'*hx2 + ones(size(z'))*x0;
whx = w'*hx2;
xx = xx(:);
whx=whx(:);

y0 = (y(1:end-1)+y(2:end))*0.5;
hy2 = diff(y)*0.5;
yy = z'*hy2 + ones(size(z'))*y0;
why = w'*hy2;
yy = yy(:);
why=why(:);

[XX,YY]=meshgrid(xx,yy);
ZZ=griddata(X',Y',fun_xy,XX,YY);



[c_n] = fun_compute2DFouruerCoeff_mex(nn,mm,XX,YY,ZZ,whx,why);



% % c_n=zeros(length(mm),length(nn));
% % 
% % for ii=1:length(mm)
% %     for jj=1:length(nn)
% %         expmt = exp(-1i*(XX*mm(ii)+YY*nn(jj)))/(4*pi*pi);
% %         c_n(ii,jj)  = (ZZ.*expmt*whx)'*why;
% %     end
% % end

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
end