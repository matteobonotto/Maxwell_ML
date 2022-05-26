function [Js]=fourier_inv_fv(JJ,nharm,ndim2,m_start,theta)

Js=zeros(ndim2,1);
for kharm=1:nharm
    for ks=1:ndim2
        exp1= cos((m_start+kharm-1)*theta(ks))+1i*...
            sin((m_start+kharm-1)*theta(ks));
        Js(ks)=Js(ks)+JJ(kharm)*exp1;
    end
end

end

