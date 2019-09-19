function noise = wgn(m, n, power)
    noise = power * randn(m,n);
end

