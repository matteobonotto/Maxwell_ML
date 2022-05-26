function [c_nm] = fun_DFT2D_Matlab_new(f,nm_harm)

mm_Matlab = -floor(size(f,1)/2):floor(size(f,1)/2);
nn_Matlab = -floor(size(f,2)/2):floor(size(f,2)/2);

qq = fftshift(fft2(f))/numel(f);

c_nm = zeros(size(nm_harm,1),1);
for ii = 1:size(nm_harm,1)
    c_nm(ii) = qq(mm_Matlab == nm_harm(ii,2), nn_Matlab == nm_harm(ii,1));
end
