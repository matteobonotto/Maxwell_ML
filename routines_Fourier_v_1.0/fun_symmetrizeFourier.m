function [c_mn_new,mm_new,nn_new]=fun_symmetrizeFourier(c_m,mm,nn)

if nn<0
    nn_new=[nn -nn];
elseif nn>0
    nn_new=[-nn nn];
end

if -mm(1)==mm(end) %ALREADY SYMMETRIC
    mm_new=mm;
    nM_new=length(mm_new);
    c_mn_new=zeros(nM_new,2);

    c_mn_new(:,1)=c_m;
    c_mn_new(:,2)=conj(flipud(c_mn_new(:,1)));
        
else
    mm_new=-max(abs(mm)):max(abs(mm));
    nM_new=length(mm_new);
    c_mn_new=zeros(nM_new,2);
    
    % SELEZIONO PARTE SIMMETRICA DELLE ARMONICHE
    imm_symm=1:find(mm==-mm(1));
    mm_symm=mm(imm_symm);
    
    i_1=find(mm_new==mm_symm(1));
    i_2=find(mm_new==mm_symm(end));
    
    c_mn_new(i_1:i_2,1)=c_m(imm_symm,1);
    
    % SELEZIONO LE RIMANENTI ARMONICHE
    imm_asymm=imm_symm(end)+1:length(mm);
    c_mn_new(i_2+1:end,1)=c_m(imm_asymm,1);
    
    % GENERO LE ARMONICHE MANCANTI PER N SIMMETRICO
    c_mn_new(:,2)=flipud(conj(c_mn_new(:,1)));
end

end