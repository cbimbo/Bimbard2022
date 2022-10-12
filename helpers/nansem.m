function sem = nansem(x)
    %%% Computes the standard error of the mean.

    x = x(~isnan(x));
    sem = std(x)/sqrt(numel(x));