function z = ztransform(c,type)

    if ~exist('type','var')
        type = 'none';
    end

    switch type
        case 'none'
            % simple Z-transform
            z = atanh(c);
        case 'mean'
            z = tanh(nanmean(atanh(c)));
        case 'std'
            z = tanh(nanstd(atanh(c)));
        case 'sem'
            z = tanh(nanstd(atanh(c)))/sqrt(sum(~isnan(c)));
    end