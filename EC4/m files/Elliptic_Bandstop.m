clear all;
%% Unnormalized digital filter specifications (in kHz)
fp1 = 61.1;
fp2 = 89.1;
fs1 = 65.1;
fs2 = 85.1;
samp_freq = 260;

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
Omega_Ls = 1.4029;
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
analog_bsf(s) = analog_lpf((B*s)/((s*s)+(Omega_o*Omega_o)));
discrete_bsf(z) = analog_bsf((z-1)/(z+1));
%% Coefficients of Analog LPF
[nums, dens] = numden(analog_lpf(s));
num_lpf = sym2poly(expand(nums));
den_lpf = sym2poly(expand(dens));
num_lpf = num_lpf./den_lpf(1);
den_lpf = den_lpf./den_lpf(1);
[H_lpf,f_lpf] = freqs(num_lpf, den_lpf, 1000);
pos_stop = 0;
pos_start1 = 0;
for i=1:length(f_lpf)
    if(f_lpf(i)>=1.862)
        pos_stop = i;
        break;
    end
end
for i=1:length(f_lpf)
    if(f_lpf(i)>=0.9912)
        pos_start1 = i;
        break;
    end
end
for i=1:length(f_lpf)
    if(f_lpf(i)>=0.6183)
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
[nums_b, dens_b] = numden(analog_bsf(s));
num_bsf = sym2poly(nums_b);
den_bsf = sym2poly(dens_b);
num_bsf = num_bsf./den_bsf(1);
den_bsf = den_bsf./den_bsf(1);
[H_bsf,f_bsf] = freqs(num_bsf, den_bsf, 100000);
critical_points_bsf = zeros(6);
check_points_bsf = [0.9099, 1.073, 1.2966, 1.574, 1.8561, 2.262];
for i=1:6
    for l=1:length(f_bsf)
        if(f_bsf(l)>=check_points_bsf(i))
            critical_points_bsf(i) = l;
            break;
        end
    end
end
figure
hold on;
plot(f_bsf, abs(H_bsf)); 
title('H_{analog, BSF}(s) on normalized frequency axis'); 
xlim([0.1,3.5]); ylim([0, 1.2]); xlabel('s'); ylabel('|H_{analog,BSF}(s)|'); 
for i=1:6
    plot(f_bsf(critical_points_bsf(i)), abs(H_bsf(critical_points_bsf(i))),'r*');
end
grid on;
%% Discrete Filter Coefficients

[nums_b2, dens_b2] = numden(discrete_bsf(z));
num_bsf2 = sym2poly(nums_b2);
den_bsf2 = sym2poly(dens_b2);
num_bsf2 = num_bsf2./den_bsf2(1);
den_bsf2 = den_bsf2./den_bsf2(1);
[H_bsf2,f_bsf2] = freqz(num_bsf2, den_bsf2, 1024*1024, 260e3);
critical_points_bsf = zeros(4);
check_points_bsf = [61.1e3, 67.84e3, 83.27e3, 89.1e3];
for i=1:length(critical_points_bsf)
    for l=1:length(f_bsf2)
        if(f_bsf2(l)>=check_points_bsf(i))
            critical_points_bsf(i) = l;
            break;
        end
    end
end
figure
hold on;
plot(f_bsf2, abs(H_bsf2)); 
title('H_{discrete, BSF}(z) on un-normalized frequency axis'); 
ylim([0, 1.2]); xlabel('Frequency in Hz'); ylabel('|H_{discrete,BSF}(z)|'); 
for i=1:length(critical_points_bsf)
    plot(f_bsf2(critical_points_bsf(i)), abs(H_bsf2(critical_points_bsf(i))),'r*');
end
grid on;
%%
fvtool(num_bsf2, den_bsf2);
