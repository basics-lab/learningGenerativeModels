function out = lognormcdf (x)
    out = double(py.scipy.stats.norm().logcdf(x));
end