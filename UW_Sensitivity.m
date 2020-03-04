function new_spec = UW_Sensitivity(f, Xss, sens, source, mic) 

% Requires inputs of frequency array and single-sided spectrum from
% my_autospec function. Also requires the sensitivity value put in
% AFR.Input also the transmitter and receiver device numbers (e.g. 4034)
% Outputs the correct pressure values from recording. 

% TODO: Add ability to change what mics were used for reciever and source

%% Load Sensitivity Data

hydrophone_4034 = load('4034 Hydrophone sensitivity.csv');
hydrophone_4038 = load('4038 Hydrophone sensitivity.csv');
projector_4034 = load('4034 Projector sensitivity.csv');
projector_4038 = load('4038 Projector sensitivity.csv');


%% Adjust data to sensitivity
if mic == 4034
    f_mic = 1000.*hydrophone_4034(:,1);
    dB_mic = hydrophone_4034(:,2);
end

if mic == 4038
    f_mic = 1000.*hydrophone_4038(:,1);
    dB_mic = hydrophone_4038(:,2);
end

if source == 4038
    f_source = 1000.*projector_4038(:,1);
    dB_source = projector_4038(:,2);
end 


% now interpolate sensitivity data to fit data size of measurement'f1'
    dB_interp = interp1(f_mic,dB_mic,f, 'spline','extrap');

    new_spec = Xss*sens./(1e9.*10.^(dB_interp./20));
