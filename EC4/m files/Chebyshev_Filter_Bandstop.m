clear all;
%% Specifications of Butterworth Filter
delta = 0.15;
D1 = (1/((1-delta)^2))-1;
epsilon = sqrt(D1); %Cutoff frequency
N = 4; %Order of filter

%% Poles of Butterworth Filter 
p1 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p2 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p3 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);
p4 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);

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

%% Creating the Transfer function for the Analog Low Pass Filter
k = (p1*p2*p3*p4)/sqrt(1+D1);
[numerator, denominator] = zp2tf([], [p1 p2 p3 p4], k); % Transfer Function is multiplied with a constant to make DC gain=1 at s=0

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
[H_lpf,f_lpf] = freqs(num_lpf, den_lpf, 10000);
pos_stop = 0;
pos_start = 0;
for i=1:length(f_lpf)
    if(f_lpf(i)>=1.4029)
        pos_stop = i;
        break;
    end
end
for i=1:length(f_lpf)
    if(f_lpf(i)>=0.9996)
        pos_start = i;
        break;
    end
end
figure
hold on;
plot(f_lpf, abs(H_lpf)); 
title('H_{analog, LPF}(s_L) on normalized frequency axis'); 
xlim([0.1,3.5]); ylim([0, 1.2]); xlabel('s_L'); ylabel('|H_{analog,LPF}(s_L)|'); 
plot(f_lpf(pos_stop),abs(H_lpf(pos_stop)),'r*');
plot(f_lpf(pos_start),abs(H_lpf(pos_start)),'r*');
grid on;

%% Coefficients of Analog BSF
[nums_b, dens_b] = numden(analog_bsf(s));
num_bsf = sym2poly(nums_b);
den_bsf = sym2poly(dens_b);
num_bsf = num_bsf./den_bsf(1);
den_bsf = den_bsf./den_bsf(1);
[H_bsf,f_bsf] = freqs(num_bsf, den_bsf, 100000);
critical_points_bsf = zeros(6);
check_points_bsf = [0.9099, 1.0025, 1.2966, 1.6586, 2.135, 1.8561];
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
check_points_bsf = [61.1e3, 65.1e3, 85.1e3, 89.1e3];
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
