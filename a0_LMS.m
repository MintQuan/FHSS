function [h_es] = a0_LMS(sigma_n, h, L)
    
    % ========== Chu?i hu?n luy?n
    N = 10^5;
    xt = randsrc(1,N,[-1 1]);
    y_temp = conv(xt, h);
    n = sqrt(sigma_n /2) * randn(1, N + L - 1) + 1i *sqrt(sigma_n /2) * randn(1, N + L - 1);
    yt = y_temp + n;
    
    % ======== Thu?t toán LMS
    muy = 10^-2;
    xt = [zeros(1,L-1) xt];
    w = [1  zeros(1, L - 1)]; % Kh?i t?o
    for n = 1:N
        U = xt(n : n + (L - 1));
        U = flip(U);
        e = yt(n) - U * w.';
        w = w + muy * conj(U) *e;
    end
    h_es = w; % Kênh ??c l??ng
    
end
