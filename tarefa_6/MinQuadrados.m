% estima parametros pelo calculo de minimos quadrados
function [yHat, theta] = MinQuadrados(y, u, order)
    yPsi = y(order + 1:end);    
    PSI = [];
    for i = 0:order-1
        PSI = [PSI, y(order-i:end-i-1) u(order-i:end-i-1)];
    end
    
    theta = (PSI' * PSI) \ PSI' * yPsi;
    yHat = PSI*theta;
end