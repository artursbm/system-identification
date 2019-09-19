function deriv = food_chain_ode_sim(t, odes)
%   odes(1) = x, odes(2) = y, odes(3) = z
    a=0.311;
    b=0.518;
    c=1.036;
    d=0.311;
    e=0.161;
    f=4.599;
    g=2.469;
    h=0.322;
    
    deriv1 = (odes(1)*(1-odes(1))) - ((odes(1)*odes(2))/(odes(1)+a));
    deriv2 = (-b*odes(2)) + ((c*odes(1)*odes(2))/(odes(1)+d)) - ((odes(2)*odes(3))/(odes(2)+e));
    deriv3 = (f*(odes(3)^2)) - ((g*(odes(3)^2))/(odes(2)+h));
    
    deriv = [deriv1; deriv2; deriv3];
end

