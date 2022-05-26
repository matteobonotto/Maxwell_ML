function f = fun_fourier_inv_2D(c_nm,nn,mm,PHI,THETA)

f = zeros(size(THETA));
for ii = 1:numel(mm)
    for jj = 1:numel(nn)
% %         f = f + c_nm(ii,jj)*exp(1i*nn(jj)*PHI + 1i*mm(ii)*THETA);
        f = f + c_nm(jj,ii)*exp(1i*nn(jj)*PHI + 1i*mm(ii)*THETA);
    end
end

if any(imag(f) > 1e-10)
    warning('inaginary part of inverse Fourier transform not negligible')
    f = 2*real(f);
else
    f = real(f);
end


% % f = ifftshift(f);
