function [Fp, overestimate] = computeISIViolations(st,tauR,tauC)
    %%% From Hill et al., 2011

    N = numel(st);
    T = st(end) - st(1);

    a = 2 * (tauR - tauC) * N^2 / T;
    r = sum(diff(st) <= tauR);

    if r == 0
        Fp = 0;
        overestimate = 0;
    else
        rts = roots([-1, 1, -r / a]);
        Fp = min(rts);
        overestimate = 0;
        if ~isreal(Fp) %function returns imaginary number if r is too high: overestimate number.
            overestimate = 1;
            if r < N %to not get a negative weird number or a 0 denominator
                Fp = r / (2 * (tauR - tauC) * (N - r));
            else
                Fp = 1;
            end
        end
    end

end