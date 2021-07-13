%UWIntensityCalc
%This script gives an implementation of TRAD_func in finding acoustic
%vector intensity in underwater applciations. This code also can use the
%PAGE_func method to find intensity, though this may have issues in the
%reverberant environment of the UW Lab tank.

%You can choose whether to use data or a simulated field by executing the
%appropriate segments of code

%% Constants (Run this first using CTRL+ENTER)
c = 1478; % ~1500 for water
rho = 1023; % value for water
Iref = 6.61e-19; % Iref = Pref^2/Z0 = Pref^2/(c*rho)
Pref = 1e-6; % Pref = 1e-6 for water

a = .34; % distance between 2 microphones

l = 3; % length of recording
fs = 150000;
ns = l*fs;
df = fs/ns;
fss = 0:df:(fs/2-df);

%% Using data
% To use data, run this segment of code
% For a simulated field, run the next segment of code instead

w = hann(ns)';
W = mean(w.*conj(w)); %Used to scale the ASD for energy conservation

%An example file path--use your own here
path = 'W:\uw-measurements-tank\2021-02-17\2021-02-17_scan8';
IDnum = 2;

% Load in data, make sure binfileload.m is downloaded from byuarg library.
% This is also a good place to account for the amplification factor from 
% the NEXUS conditioning amplifier (in this case divide by 10 to account 
% for the 10 mV/Pa amplification).
x1 = binfileload(path,'ID',IDnum,1)/.10; 
x2 = binfileload(path,'ID',IDnum,2)/.10; 

% Make sure computeBlockFFt.m is downloaded from byuarg library.
[Xss1, numblocks] = computeBlockFFt(x1,ns,w,W);
[Xss2, numblocks] = computeBlockFFt(x2,ns,w,W);

Xss = cat(1, reshape(Xss1, [1 size(Xss1,1) size(Xss1,2)]),...
    reshape(Xss2, [1 size(Xss1,1) size(Xss1,2)]));

%% Calculate intensities

% This code also has capabilities to use PAGE method of finding intensity.

% PAGE = PAGE_func(fss,Xss,probe_config,rho,c,'slope',8);
TRAD = TRAD_func(fss,Xss,probe_config,rho,c);

% active intensities in the y direction
% I_PAGE_y = 10*log10(abs(PAGE.I(2,:))/Iref);
I_TRAD_y = 10*log10(abs(TRAD.I(2,:))/Iref);

%reactive
% Q_PAGE_y = 10*log10(abs(PAGE.Q(2,:))/Iref);
Q_TRAD_y = 10*log10(abs(TRAD.Q(2,:))/Iref);

TRAD_mag = 10*log10(sqrt(TRAD.Q(2,:).^2 + TRAD.I(2,:).^2)/Iref);

%% Plot
myFigureDefaults(1.8, 'sm');

figure
plot(fss,I_TRAD_y,'color',[0 0.5 0])
hold on
plot(fss,Q_TRAD_y,'r',fss,TRAD_mag,'b')
title('Y-mount 41 cm From Source')
legend('Active Intensity','Reactive Intensity','Absolute Magnitude','Location','northeast')
ylabel('I (dB re. 1 \muPa)') 
xlabel('Frequency (Hz)')
ylim([10 120])
xlim([0 3000])
pbaspect([3 1 1])
hold off

% This code can be used to compare the absolute intensity values from
% various measurement positions.

% figure
% plot(fss,TRAD_mag2,'r')
% hold on
% plot(fss,TRAD_mag1,'color',[0 0.5 0])
% plot(fss,TRAD_mag,'b')
% title('Absolute Magnitudes of Intensity')
% legend('41 cm','59 cm','77 cm','Location','northeast')
% ylabel('I (dB re. 1 \muPa)') 
% xlabel('Frequency (Hz)')
% ylim([10 120])
% xlim([0 3000])
% pbaspect([3 1 1])
% hold off

%% plotting waveforms

far = x1;
time = 0:1/fs:l-1/fs;

figure()
plot(time,(near-mean(near))*.001,'color',[0 0.5 0])
hold on 
plot(time,(mid-mean(mid))*.001,'r',time,(far-mean(far))*.001,'b')
title('Recorded Signals')
legend('41 cm','59 cm','77 cm','Location','northeast')
xlim([.5 2.5])
ylabel('Pressure (Pa)') 
xlabel('Time (s)')
hold off