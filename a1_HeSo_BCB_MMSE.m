function [c] = a1_HeSo_BCB_MMSE(N, L, h, P, sigma_n)
% ======== Ma tr?n q
q = [zeros(1, N - 1)  1  zeros(1, L - 1)];
 
% ======== Ma tr?n H
H = [];
for i = 1:N
    Temp = [zeros(1, i - 1)   h   zeros(1, N - i)];
    H = [H ; Temp];
end
 
% ======== ?áp ?ng b? cân b?ng
c = P * q * H' * (P * H * H' + sigma_n * eye(N))^-1;
end
