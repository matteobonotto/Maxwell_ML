function  [c_n] = fun_compute2DFouruerCoeff(nn,mm,XX,YY,ZZ,whx,why) %#codegen

c_n=zeros(length(nn),length(mm))+0i;

% % f_xy = @(x,y) peaks(x,y); toll_value = 1e-10;

parfor ii=1:length(nn)
    
    c_n_ii = zeros(1,length(mm))+0i;
    for jj=1:length(mm)
        expmt = exp(-1i*(XX*nn(ii)+YY*mm(jj)))/(4*pi*pi);
        
        tmp = (ZZ.*expmt*whx).'*why;
        
        % %         g_ij = @(x,y) exp(-1i*nn(ii)*x -1i*mm(jj)*y);
        % %         f_ij = @(x,y) f_xy(x,y) .* g_ij(x,y);
        % %         temp_int = integral2(f_ij,-pi,pi,-pi,pi,'AbsTol',toll_value,'RelTol',toll_value);
        % %         temp_int/4/pi/pi
        
        c_n_ii(1,jj)  = tmp(1,1);
    end
    
    c_n(ii,:) = c_n_ii;
end