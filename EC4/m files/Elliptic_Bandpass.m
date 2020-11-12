clear all;
%% Unnormalized digital filter specifications (in kHz)
fp1 = 77.9e3;
fp2 = 97.9e3;
fs1 = 73.9e3;
fs2 = 101.9e3;
samp_freq = 330e3;

%% Normalized Analog Filter Specifications
omega_p1 = tan((fp1*pi)/samp_freq);
omega_p2 = tan((fp2*pi)/samp_freq);
omega_s1 = tan((fs1*pi)/samp_freq);
omega_s2 = tan((fs2*pi)/samp_freq);

%% Paramters for BandStop to LPF Transformation
Omega_o = sqrt(omega_p1*omega_p2);
B = (omega_p2-omega_p1);

%% Specifications of Low Pass Filter
delta_1 = 0.15;
delta_2 = 0.15;
D_1 = ((1/(1-delta_1))^2)-1;
D_2 = ((1/delta_2)^2)-1;
k_1 = sqrt(D_1/D_2);
Omega_Ls = 1.4057;
Omega_Lp = 1;
k = Omega_Lp/Omega_Ls;
%% Calculating the complete elliptic integral(and its complement) of first kind for k and k_1 using ellipk function
[K,Kc] = ellipk(k);
[K1,K1c] = ellipk(k_1);
%% Calcultating the minimum degree of the filter and recalculating k and k_1 for exact(and more stringent) passband and stopband characteristics
N_min = ceil((K1c*K)/(K1*Kc));
k = ellipdeg(N_min,k_1);
Omega_Ls = Omega_Lp/k; % new Omega_Ls(more stringent)
L = floor(N_min/2);
r = (N_min-(2*L));
i = (1:L)';
u = (2*i-1)/N_min;
%% Finding zeroes and poles of the LPF
zeta = cde(u,k);
zeroes_lpf = (1j)./(k*zeta);
v0 = (-1j)*asne(1j/sqrt(D_1),k_1)/N_min;
poles_lpf = 1j*cde(u-1j*v0,k);
pole_0 = 1j*sne(1j*v0,k);
%% Writing Numerator and Denominator
Constant_coeff = 1;
for i=1:L
    Constant_coeff = Constant_coeff*((abs(poles_lpf(i)/abs(zeroes_lpf(i))))^2);
end
if r==1
    Constant_coeff = Constant_coeff*(-pole_0);
else
    Constant_coeff = Constant_coeff/sqrt(1+D_1);
end
zeroes_lpf = [zeroes_lpf,conj(zeroes_lpf)]';
poles_lpf = [poles_lpf,conj(poles_lpf),pole_0]';
[numerator, denominator] = zp2tf(zeroes_lpf, poles_lpf, Constant_coeff);
%% Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(numerator,s)/poly2sym(denominator,s);
analog_bpf(s) = analog_lpf(((s*s)+(Omega_o*Omega_o))/(B*s));
discrete_bpf(z) = analog_bpf((z-1)/(z+1));
%% Coefficients of Analog LPF
[nums, dens] = numden(analog_lpf(s));
num_lpf = sym2poly(expand(nums));
den_lpf = sym2poly(expand(dens));
num_lpf = num_lpf./den_lpf(1);
den_lpf = den_lpf./den_lpf(1);
[H_lpf,f_lpf] = freqs(num_lpf, den_lpf, 10000);
pos_stop = 0;
pos_start1 = 0;
pos_start2 = 0;
for i=1:length(f_lpf)
    if(f_lpf(i)>=1.875)
        pos_stop = i;
        break;
    end
end
for i=1:length(f_lpf)
    if(f_lpf(i)>=0.999)
        pos_start1 = i;
        break;
    end
end
for i=1:length(f_lpf)
    if(f_lpf(i)>=0.6095)
        pos_start2 = i;
        break;
    end
end
figure
hold on;
plot(f_lpf, abs(H_lpf)); 
title('H_{analog, LPF}(s_L) on normalized frequency axis'); 
xlim([0.1,3.5]); ylim([0, 1.2]); xlabel('s_L'); ylabel('|H_{analog,LPF}(s_L)|'); 
plot(f_lpf(pos_stop),abs(H_lpf(pos_stop)),'r*');
plot(f_lpf(pos_start1),abs(H_lpf(pos_start1)),'r*');
plot(f_lpf(pos_start2),abs(H_lpf(pos_start2)),'r*');
grid on;
%% Coefficients of Analog BSF
[nums_b, dens_b] = numden(analog_bpf(s));
num_bpf = sym2poly(nums_b);
den_bpf = sym2poly(dens_b);
num_bpf = num_bpf./den_bpf(1);
den_bpf = den_bpf./den_bpf(1);
[H_bpf,f_bpf] = freqs(num_bpf, den_bpf, 10000);
critical_points_bpf = zeros(5);
check_points_bpf = [0.7765, 0.91599, 0.9793, 1.342, 1.583];
for i=1:5
    for l=1:length(f_bpf)
        if(f_bpf(l)>=check_points_bpf(i))
            critical_points_bpf(i) = l;
            break;
        end
    end
end
figure
hold on;
plot(f_bpf, abs(H_bpf)); 
title('H_{analog, BPF}(s) on normalized frequency axis'); 
xlim([0.1,3.5]); ylim([0, 1.2]); xlabel('s'); ylabel('|H_{analog,BPF}(s)|'); 
for i=1:5
    plot(f_bpf(critical_points_bpf(i)), abs(H_bpf(critical_points_bpf(i))),'r*');
end
grid on;
%% Discrete Filter Coefficients

[nums_b2, dens_b2] = numden(discrete_bpf(z));
num_bpf2 = sym2poly(nums_b2);
den_bpf2 = sym2poly(dens_b2);
num_bpf2 = num_bpf2./den_bpf2(1);
den_bpf2 = den_bpf2./den_bpf2(1);
[H_bpf2,f_bpf2] = freqz(num_bpf2, den_bpf2, 1024*1024, samp_freq);
critical_points_bpf = zeros(4);
check_points_bpf = [69.68e3, 97.9e3, 77.9e3, 105.8e3];
for i=1:length(critical_points_bpf)
    for l=1:length(f_bpf2)
        if(f_bpf2(l)>=check_points_bpf(i))
            critical_points_bpf(i) = l;
            break;
        end
    end
end
figure
hold on;
plot(f_bpf2, abs(H_bpf2)); 
title('H_{discrete, BSF}(z) on un-normalized frequency axis'); 
ylim([0, 1.2]); xlabel('Frequency in Hz'); ylabel('|H_{discrete,BSF}(z)|'); 
for i=1:length(critical_points_bpf)
    plot(f_bpf2(critical_points_bpf(i)), abs(H_bpf2(critical_points_bpf(i))),'r*');
end
grid on;
%%
fvtool(num_bpf2, den_bpf2);
%%
plot(real(poles_lpf),imag(poles_lpf),'o');
