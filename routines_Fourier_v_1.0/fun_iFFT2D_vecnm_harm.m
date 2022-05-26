function f = fun_iFFT2D_vecnm_harm(c_nm,nm_harm,PHI,THETA)


nn = nm_harm(:,1);
mm = nm_harm(:,2);

f = zeros(size(THETA));
for ii = 1:numel(mm)
    f = f + c_nm(ii)*exp(1i*nn(ii)*PHI + 1i*mm(ii)*THETA);
end

if any(imag(f(:)) > 1e-10)
    % %     warning('inaginary part of inverse Fourier transform not negligible')
    f = 2*real(f);
else
    f = real(f);
end












