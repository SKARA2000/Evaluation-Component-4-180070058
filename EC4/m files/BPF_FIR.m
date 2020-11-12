clear all;
%%
f_samp = 330e3;
delta = 0.15; %stopband tolerance
A = -20*log10(delta);

%% Shape parameter
alpha=0;
if(A<21)
    alpha=0;
elseif((A>=21) && (A<=50))
    alpha = 0.5842*((A-21)^0.4)+0.07886*(A-21);
else
    alpha = 0.1102*(A-8.7);
end
%% Length of FIR filter
delta_omega = 0.0242*pi;
N = ceil((A-8)/(2*2.285*delta_omega));

%% Kaiser Window formation
N_corrected = (N+9); %N is corrected due to poor lower bound derived previously
n = (2*(N_corrected)+1);
beta = alpha/N;
kaiser_window = kaiser(n,beta);

%% FIR BPF Filter formation
passband_1 = 0.46*pi;
passband_2 = 0.60545*pi;
bpf = lpf_FIR(passband_2,n) - lpf_FIR(passband_1,n);
bpf_FIR = bpf.*kaiser_window';
fvtool(bpf_FIR);
%% Plotting FIR filter on unnormalized frequencies
[H_FIR, f_FIR] = freqz(bpf_FIR,1,1024,f_samp);
critical_points_bpf = zeros(6);
check_points_bpf = [77.9e3, 97.9e3, 73.9e3, 101.9e3, 95.07e3, 80.57e3];
for i=1:length(critical_points_bpf)
    for l=1:length(f_FIR)
        if(f_FIR(l)>=check_points_bpf(i))
            critical_points_bpf(i) = l;
            break;
        end
    end
end
figure
hold on;
plot(f_FIR, abs(H_FIR)); 
title('H_{discrete, BPF}(z) on un-normalized frequency axis'); 
ylim([0, 1.2]); xlabel('Frequency in Hz'); ylabel('|H_{discrete,BPF}(z)|'); 
for i=1:length(critical_points_bpf)
    plot(f_FIR(critical_points_bpf(i)), abs(H_FIR(critical_points_bpf(i))),'r*');
end
grid on;
