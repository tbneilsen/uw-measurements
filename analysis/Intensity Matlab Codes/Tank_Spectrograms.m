%This code creates spectrograms from data taken in BYU UW Lab Tank. As such
%the acoustic parameters used in the calculations are correct for in-water
%measurements. 

clear all; close all;

%% Recording Parameters and Load in Data

path = 'W:\uw-measurements-tank\2021-02—11\2021-02-11_scan15';

CHnum = 2; % the channel of interest from input to DAQ 
fs = 150000; % sampling frequency
length = 3; % time length of signal, including preceeding and trailing 0's

% Load in data, make sure binfileload.m is downloaded from byuarg library.
% Also, set IDnum for the ID of interest. This is also a good place to
% account for the amplification factor from the NEXUS conditioning
% amplifier (in this case divide by 10 to account for the 10 mV/Pa
% amplification).
data(1,:) = binfileload( path1, 'ID',0,CHnum)./10;
data(2,:) = binfileload( path1, 'ID',6,CHnum)./10;
data(3,:) = binfileload( path1, 'ID',12,CHnum)./10;
data(4,:) = binfileload( path1, 'ID',18,CHnum)./10;
data(5,:) = binfileload( path1, 'ID',24,CHnum)./10;

% set a time array of same length as the loaded data
time = 0:1/fs:length-1/fs;

%% Plotting Time

% Make sure specgram.m is downloaded from byuarg library. Mean is
% subtracted from data to remove any possible DC offset variation between
% measurements. Last parameter, ns, should be set so that there is
% sufficient time resolution and frequency bin size for your purposes.
% Generally ns=15000 for spectrograms
[Gxx,t,f,runOASPL] = specgram(data(1,:) - mean(data(1,:)),fs,15000); 
[Gxx2,t,f,runOASPL2] = specgram(data(2,:)-mean(data(2,:)),fs,15000);
[Gxx3,t,f,runOASPL3] = specgram(data(3,:)-mean(data(3,:)),fs,15000);
[Gxx4,t,f,runOASPL4] = specgram(data(4,:)-mean(data(4,:)),fs,15000);
[Gxx5,t,f,runOASPL5] = specgram(data(5,:)-mean(data(5,:)),fs,15000);

% See instructions in myFigureDefaults.m from byuarg library.
myFigureDefaults(1.8, 'sm');

figure()
pcolor(t,f,10*log10(Gxx5'/1e-6))
ylim([0 3000])
caxis([-20 10])
title('Near Anechoic Lining')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
shading interp
c = colorbar;
c.Label.String = 'SPL (dB re. 1 \muPa)';


figure()
pcolor(t,f,10*log10(Gxx3'/1e-6))
ylim([0 3000])
caxis([-20 10])
title('Middle')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
shading interp
c = colorbar;
c.Label.String = 'SPL (dB re. 1 \muPa)';


figure()
pcolor(t,f,10*log10(Gxx'/1e-6))
ylim([0 3000])
caxis([-20 10])
title('Near Tank Wall')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
shading interp
c = colorbar;
c.Label.String = 'SPL (dB re. 1 \muPa)';

%% Line Plots of Acoustic Energy at Specific Frequencies

myFigureDefaults(1.8, 'sm');

% These plots illustrate the amount of acoustical energy at the specified
% frequency across the whole time of the signal. Checking the specified
% index number in the frequency array f(x) will show at what frequency the
% plot is given.
figure()
plot(t,20*log10(Gxx(:,6)/1e-12),'k', t,20*log10(Gxx5(:,6)/1e-12),'r-.')
xlabel('Time (s)')
ylabel('SPL (dB re. 1 \muPa)')
grid on
title('SPL at XX Hz')
legend({'near wall', 'near anechoic'})

figure()

plot(t,20*log10(Gxx(:,11)/1e-12),'k',t,20*log10(Gxx5(:,11)/1e-12),'r-.')
xlabel('Time (s)')
ylabel('SPL (dB re. 1 \muPa)')
grid on
title('SPL at YY Hz')
legend({'near wall', 'near anechoic'})

%% Make Schematics of Experimental Set-up

figure()
plot([.2,.2],[2.36,2.14],'r*'...
    ,[.6,.6],[2.36,2.14],'g*',[1,1],[2.36,2.14],'b*')
xlim([0 1.2])
ylim([0 3.6])
patch([1.14 1.14 1.2 1.2], [1.5 3 3 1.5], [0,0,1])
legend('D','E','F','Location','southwest')
pbaspect([1 3 1])
