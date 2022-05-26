function [c_nm,NN,MM] = fun_DFT2D_Matlab(f,nn,mm)

mm_Matlab = -floor(size(f,1)/2):floor(size(f,1)/2);
nn_Matlab = -floor(size(f,2)/2):floor(size(f,2)/2);

[NN,MM] = meshgrid(nn,mm);

qq = fftshift(fft2(f))/numel(f);

c_nm = zeros(numel(mm),numel(nn));
for ii = 1:numel(mm)
    for jj = 1:numel(nn)
        c_nm(ii,jj) = qq(mm_Matlab == mm(ii), nn_Matlab == nn(jj));
    end
end
