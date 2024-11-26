clc
clear

% Generation of Sample bits
sequence=round(rand(1,20));    % Generating 20 bits
input_signal=[];                 % input signal declaration
carrier_signal=[];               % carrier signal declaration
time= 0:2*pi/119:2*pi;          % Number of samples - Creating a vector of 120 total values including end point [0:0.052:6.28]
L = 5; %số kênh đa đường
N = 20; %chiều dài của bộ cân bằng
Ns = 20;
sigma_h = 1;
sigma_n = 0.0001;
P = 10^5;
N_samples = 120;

for k = 1 :Ns
    if sequence(1,k)==0
        sig=-ones(1,120);    % -1 value for binary input 0
    else
        sig=ones(1,120);     % +1 value for binary input 1
    end
    c=cos(time);             % Carrier frequency is calculated.
    carrier_signal = [carrier_signal c];  %signal that is generated which will transmitted with input signal
    input_signal = [input_signal sig];    %Original signal to be transmitted
end

figure(1)
subplot(4,1,1);
plot(input_signal);    % Plotting input signal
axis([-100 2400 -1.5 1.5]);
title('\bf\it Original 20 bit Sequence');

% BPSK Modulation of the signal
bpsk_mod_signal=input_signal.*carrier_signal;   % Modulating the signal
subplot(4,1,2);
plot(bpsk_mod_signal);  %Plotting BPSK Modulated signal
axis([-100 2400 -1.5 1.5]);
title('\bf\it BPSK Modulated Signal');

% Preparation of 6 new carrier frequencies
%6 frequencies are randomly selected by the PN sequence generator
time1=[0:2*pi/9:2*pi];  %f1=13.33
time2=[0:2*pi/19:2*pi]; %f2 = 6.32
time3=[0:2*pi/29:2*pi]; %f3 = 4.14
time4=[0:2*pi/39:2*pi]; %f4 = 3.08
time5=[0:2*pi/59:2*pi]; %f5 = 2.03
time6=[0:2*pi/119:2*pi];%f6 = 1.01
carrier1=cos(time1);
carrier1=[carrier1 carrier1 carrier1 carrier1 carrier1 carrier1 carrier1 carrier1 carrier1 carrier1 carrier1 carrier1];
carrier2=cos(time2);
carrier2=[carrier2 carrier2 carrier2 carrier2 carrier2 carrier2];
carrier3=cos(time3);
carrier3=[carrier3 carrier3 carrier3 carrier3];
carrier4=cos(time4);
carrier4=[carrier4 carrier4 carrier4];
carrier5=cos(time5);
carrier5=[carrier5 carrier5];
carrier6=cos(time6);

% Random frequency hopps to form a spread signal
spread_signal=[];
for n=1:20
    c=randi([1 6],1,1);     %6 frequencies are randomly are selected by the PN generator
    switch(c)
        case(1)
            spread_signal=[spread_signal carrier1];
        case(2)
            spread_signal=[spread_signal carrier2];
        case(3)
            spread_signal=[spread_signal carrier3];
        case(4)
            spread_signal=[spread_signal carrier4];
        case(5)
            spread_signal=[spread_signal carrier5];
        case(6)
            spread_signal=[spread_signal carrier6];
    end
end
subplot(4,1,3)
plot([1:2400],spread_signal);
axis([-100 2400 -1.5 1.5]);
title('\bf\it Spread Signal with 6 frequencies');

% Spreading BPSK Signal
freq_hopped_sig=bpsk_mod_signal.*spread_signal;    %This is the signal which is being finally transmitted
subplot(4,1,4)
plot([1:2400],freq_hopped_sig);
axis([-100 2400 -1.5 1.5]);
title('\bf\it Frequency Hopped Spread Spectrum Signal');

% Kênh truyền
h = sqrt(sigma_h/2) * (randn(1,L) + 1i*sqrt(sigma_h/2)*randn(1,L));

%Kênh nhiễu
n = sqrt(sigma_n/2)*randn(1,length(freq_hopped_sig)+L-1) + 1i*sqrt(sigma_n/2)*randn(1,length(freq_hopped_sig)+L-1);

% Khối thu
y = conv(freq_hopped_sig,h) + n;
y = [y zeros(1,N-L)];
% Ước lượng kênh
[h_es] = a0_LMS(sigma_n, h, L);
c_MMSE = a1_HeSo_BCB_MMSE(N, L, h_es, P, sigma_n);

z = zeros(1,120*Ns);
for j = 1:120*Ns
    Y = y(:,j:j+Ns-1);
    Y = transpose(flip(Y));
    z(j) = c_MMSE*Y;
end

bpsk_demodulated = real(z)./spread_signal;
figure(2)
subplot(2,1,1)
plot([1:2400],bpsk_demodulated);
axis([-100 2400 -1.5 1.5]);
title('\bf Demodulated BPSK Signal from Wide Spread');
original_BPSK_signal=bpsk_demodulated./carrier_signal;     %FFH demodulated signal is data demodulated by means of BPSK Signal
% Phân chia tín hiệu thành các khung 120 mẫu
decoded_bits = [];  % Mảng lưu kết quả giải mã
% Phân chia tín hiệu thành các khung 120 mẫu
decoded_bits = [];  % Mảng lưu kết quả giải mã
for i = 1:Ns
    start_idx = (i-1) * N_samples + 1;
    end_idx = i * N_samples;
    bit_signal = original_BPSK_signal(start_idx:end_idx);

    % Quy tắc Majority Voting (Bỏ phiếu đa số)
    % Đếm số mẫu gần 1 và gần 0
    num_zeros = sum(bit_signal < 0);  % Đếm số mẫu gần 0
    num_ones = N_samples - num_zeros;  % Đếm số mẫu gần 1

    % Quyết định bit dựa trên số lượng mẫu gần 0 hoặc gần 1
    if num_zeros > num_ones
        decoded_bits = [decoded_bits, 0];  % Nhận bit 0
    else
        decoded_bits = [decoded_bits, 1];  % Nhận bit 1
    end
end
temp = [];
for i = 1:Ns
    if decoded_bits(1,i) == 0
        rx = -ones(1,120);
    else
        rx = ones(1,120);
    end
temp = [temp rx];
end
decoded_bits = temp;
subplot(2,1,2)
plot([1:2400],decoded_bits);
axis([-100 2400 -1.5 1.5]);
title('\bf Transmitted Original Bit Sequence');




