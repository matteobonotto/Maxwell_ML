function [Js]=fun_FFT_inv(c_n,mm,theta)

m_start = mm(1);
nharm = length(mm);
ndim2 = length(theta);

Js=zeros(ndim2,1);
for kharm=1:nharm
    for ks=1:ndim2
        exp1= cos((m_start+kharm-1)*theta(ks))+1i*...
            sin((m_start+kharm-1)*theta(ks));
        Js(ks)=Js(ks)+c_n(kharm)*exp1;
    end
end


if any(imag(Js) > 1e-10)
    warning('inaginary part of inverse Fourier transform not negligible')
    Js = 2*real(Js);
else
    Js = real(Js);
end
 


end

