function f = fun_fourier_inv_2D_new(c_nm,nm_harm,PHI,THETA)

f = zeros(size(THETA));
for ii = 1:size(nm_harm,1)
    f = f + c_nm(ii)*exp(1i*nm_harm(ii,1)*PHI + 1i*nm_harm(ii,2)*THETA);
end

if any(imag(f) > 1e-10)
    warning('inaginary part of inverse Fourier transform not negligible')
    f = 2*real(f);
else
    f = real(f);
end

% % f = ifftshift(f);
