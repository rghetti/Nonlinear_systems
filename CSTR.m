function out = CSTR(x, Ts, teta, K, xf, M, alpha, xc, u)

    out1 = x(1) + Ts/teta * (1 - x(1)) - Ts * K * x(1) * exp(-M / x(2));
    out2 = x(2) + Ts/teta * (xf - x(2)) + Ts * K * x(1) * exp(-M / x(2)) - Ts * alpha * u * (x(2) - xc);
    out = [out1; out2];
end