clear all;
%% Specifications of Butterworth Filter
OmegaC = 1.0722; %Cutoff frequency
N = 7; %Order of filter

%% Poles of Butterworth Filter 
p1 = OmegaC*cos(pi/2 + pi/14) + i*OmegaC*sin(pi/2 + pi/14);
p2 = OmegaC*cos(pi/2 + pi/14) - i*OmegaC*sin(pi/2 + pi/14);
p3 = OmegaC*cos(pi/2 + pi/14+pi/7) + i*OmegaC*sin(pi/2 + pi/14+pi/7);
p4 = OmegaC*cos(pi/2 + pi/14+pi/7) - i*OmegaC*sin(pi/2 + pi/14+pi/7);
p5 = OmegaC*cos(pi/2 + pi/14+2*pi/7) + i*OmegaC*sin(pi/2 + pi/14+2*pi/7);
p6 = OmegaC*cos(pi/2 + pi/14+2*pi/7) - i*OmegaC*sin(pi/2 + pi/14+2*pi/7);
p7 = -OmegaC;

%% Unnormalized digital filter specifications (in kHz)
fp1 = 77.9;
fp2 = 97.9;
fs1 = 73.9;
fs2 = 101.9;
samp_freq = 330;

%% Normalized Analog Filter Specifications
omega_p1 = tan((fp1*pi)/samp_freq);
omega_p2 = tan((fp2*pi)/samp_freq);
omega_s1 = tan((fs1*pi)/samp_freq);
omega_s2 = tan((fs2*pi)/samp_freq);

%% Paramters for BandStop to LPF Transformation
Omega_o = sqrt(omega_p1*omega_p2);
B = (omega_p2-omega_p1);

%% Creating the Transfer function for the Analog Low Pass Filter
[numerator, denominator] = zp2tf([], [p1 p2 p3 p4 p5 p6 p7], OmegaC^N); % Transfer Function is multiplied with a constant to make DC gain=1 at s=0

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
pos_start = 0;
for i=1:length(f_lpf)
    if(f_lpf(i)>=1.4057)
        pos_stop = i;
        break;
    end
end
for i=1:length(f_lpf)
    if(f_lpf(i)>=1)
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

%% Coefficients of Analog LPF
[nums_b, dens_b] = numden(analog_bpf(s));
num_bpf = sym2poly(nums_b);
den_bpf = sym2poly(dens_b);
num_bpf = num_bpf./den_bpf(1);
den_bpf = den_bpf./den_bpf(1);
[H_bpf,f_bpf] = freqs(num_bpf, den_bpf, 100000);
critical_points_bpf = zeros(5);
check_points_bpf = [0.8484, 0.9159, 1.1106, 1.342, 1.4959];
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
[H_bpf2,f_bpf2] = freqz(num_bpf2, den_bpf2, 1024*1024, 330e3);
critical_points_bpf = zeros(4);
check_points_bpf = [77.9e3, 97.9e3, 73.9e3, 101.9e3];
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
title('H_{discrete, BPF}(z) on un-normalized frequency axis'); 
ylim([0, 1.2]); xlabel('Frequency in Hz'); ylabel('|H_{discrete,BPF}(z)|'); 
for i=1:length(critical_points_bpf)
    plot(f_bpf2(critical_points_bpf(i)), abs(H_bpf2(critical_points_bpf(i))),'r*');
end
grid on;
%%
fvtool(num_bpf2, den_bpf2);