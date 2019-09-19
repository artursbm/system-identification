function y = model_2(x, k)
    
    y =  3.38040*x(k-1) - 4.30812*x(k-2) + 2.56162*x(k-3) - 1.06161*x(k-4);
    y = y - 1.21955*(x(k-5)^2) + 2.56978*x(k-1)*x(k-5)*x(k-6);
    y = y - 3.26196*x(k-3)*x(k-4)*x(k-6) + 0.48632*x(k-5);
    y = y + 2.53047*(x(k-4)^2)*x(k-5) + 0.80920*x(k-4)*x(k-7);
    y = y - 0.00455223*(x(k-1)^2)*1 + 1.47483*x(k-3)*x(k-6);
    y = y - 0.23716*(x(k-5)^2)*x(k-6) - 0.74444*x(k-1)*x(k-7);
    y = y - 0.45312*(x(k-6)^2) + 0.50283*(x(k-2)^2)*x(k-3);
    y = y - 2.02429*x(k-1)*x(k-4)*x(k-5);

end
