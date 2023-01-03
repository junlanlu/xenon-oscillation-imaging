function data = inverseboxcox(lambda, transdata, scale_factor)
    data = (transdata*lambda+1).^(1./lambda) - scale_factor;
end