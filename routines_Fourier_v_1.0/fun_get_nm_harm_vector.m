function nm_harm = fun_get_nm_harm_vector(nn,mm)

[NN,MM] = meshgrid(nn,mm);

nm_harm = [NN(:) MM(:)];


