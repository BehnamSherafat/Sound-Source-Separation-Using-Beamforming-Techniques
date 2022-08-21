%% 
tic;
clear;
close all;
clc;

%% Read data

% [s1, fs1] = audioread('BIGD9GBulldozer.wav'); % fmax = 7200 Hz
% [s2, fs2] = audioread('Jackhammer.wav');      % fmax = 132000 Hz
% [s3, fs3] = audioread('Excavator.wav');        % fmax = 7000 Hz
% 
% s1 = s1(:,1);
% s2 = s2(:,1);
% s3 = s3(:,1);
% 
% fs = fs1;
% duration = 30;
% t = (0:1/fs:duration-1/fs);
% s1 = s1(1:fs*duration);
% s2 = s2(1:fs*duration);
% s3 = s3(1:fs*duration);

%% Denoising: method = 'mcra2'
% s1(s1==0)=10^-4;
% s2(s2==0)=10^-4;
% s3(s3==0)=10^-4;

% specsub_ns('BIGD9GBulldozer.wav', 'mcra2', 'BIGD9GBulldozer_den.wav')
% specsub_ns('Jackhammer.wav', 'mcra2', 'Jackhammer_den.wav')
% specsub_ns('Excavator.wav', 'mcra2', 'Excavator_den.wav')

% Read Modified Audio
[s1, fs1] = audioread('Loader.mp3');
[s2, fs2] = audioread('Truck.mp3');

s1 = s1(:,1);
s2 = s2(:,1);
% s3 = s3(:,1);

fs = fs1;
duration = 8;
t = (0:1/fs:duration-1/fs);
s1 = s1(1:fs*duration);
s2 = s2(1:fs*duration);
s3 = s1 + s2;
audiowrite('Loader1.wav',s1,fs2);
audiowrite('Truck1.wav',s2,fs2);
audiowrite('Mixed.wav',s3,fs2);

% s3 = s3(1:fs*duration);
%% Plot data
% figure(1);
% subplot(3,1,1);
% plot(t,s1)
% title('Bulldozer', 'FontSize', 28)
% xlabel('Time (s)', 'FontSize', 26)
% ylabel('Amplitude', 'FontSize', 26)
% set(gca,'fontsize',20)
% axis tight
% 
% subplot(3,1,2);
% plot(t,s2)
% title('Jackhammer', 'FontSize', 28)
% xlabel('Time (s)', 'FontSize', 26)
% ylabel('Amplitude', 'FontSize', 26)
% set(gca,'fontsize',20)
% axis tight
% 
% subplot(3,1,3);
% plot(t,s3)
% title('Excavator', 'FontSize', 28)
% xlabel('Time (s)', 'FontSize', 26)
% ylabel('Amplitude', 'FontSize', 26)
% set(gca,'fontsize',20)
% axis tight

% %% Plot data
% figure(2);
% subplot(3,1,1);
% y = fft(s1);                               % Compute DFT of x
% L = size(s1, 1);
% P2 = abs(y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = fs*(0:(L/2))/L;
% plot(f,P1) 
% title('Bulldozer', 'FontSize', 28)
% xlabel('f (Hz)', 'FontSize', 26)
% ylabel('|P1(f)|', 'FontSize', 26)
% set(gca,'fontsize',20)
% axis tight
% 
% subplot(3,1,2);
% y = fft(s2);                               % Compute DFT of x
% L = size(s2, 1);
% P2 = abs(y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = fs*(0:(L/2))/L;
% plot(f,P1)
% title('Jackhammer', 'FontSize', 28)
% xlabel('f (Hz)', 'FontSize', 26)
% ylabel('|P1(f)|', 'FontSize', 26)
% set(gca,'fontsize',20)
% axis tight
% 
% subplot(3,1,3);
% y = fft(s3);                               % Compute DFT of x
% L = size(s3, 1);
% P2 = abs(y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = fs*(0:(L/2))/L;
% plot(f,P1) 
% title('Excavator', 'FontSize', 28)
% xlabel('f (Hz)', 'FontSize', 26)
% ylabel('|P1(f)|', 'FontSize', 26)
% set(gca,'fontsize',20)
% axis tight

%% Inputs
incidentAngle_jackhammer = [-45;0];
incidentAngle_bulldozer = [90;0];
incidentAngle_vibrator = [180;0];

temperatureC = 20.0; % ????????????????????????????????????????
c = 331.4*sqrt(1.0+(temperatureC/273));
Nele = 6; % ???????????????????????????????????????????????????
oFeq = 7e8; % ?????????????????????????????????????????????????
radius = 0.045; % ?????????????????????????????????????????????
min_freqrange = 20; % ?????????????????????????????????????????
max_freqrange = fs/2; % ???????????????????????????????????????
desiredAngle = incidentAngle_bulldozer; % ?????????????????????
numsubbands = 59; % ???????????????????????????????????????????

microphone = phased.OmnidirectionalMicrophoneElement('FrequencyRange',[min_freqrange max_freqrange]);
array = phased.UCA('Element',microphone,'NumElements',Nele,'Radius',radius, 'ArrayNormal', 'z');
collector = phased.WidebandCollector('Sensor',array,'PropagationSpeed',c,'SampleRate',fs,'NumSubbands',numsubbands,'ModulatedInput', false);

%% preallocate
NSampPerFrame = 88200; % ??????????????????????????????????????
NTSample = duration*fs;
sigArray = zeros(NTSample,Nele);

sound_bulldozer = zeros(NTSample,1);
sound_jackhammer = zeros(NTSample,1);
sound_vibrator = zeros(NTSample,1);

% set up audio device writer
audioWriter = audioDeviceWriter('SampleRate',fs, ...
        'SupportVariableSizeInput', true);
isAudioSupported = (length(getAudioDevices(audioWriter))>1);

jackhammerFileReader = dsp.AudioFileReader('Loader.wav',...
    'SamplesPerFrame',NSampPerFrame);
bulldozerFileReader = dsp.AudioFileReader('Truck.wav',...
    'SamplesPerFrame',NSampPerFrame);
vibratorFileReader = dsp.AudioFileReader('Excavator.wav',...
    'SamplesPerFrame',NSampPerFrame);

%% simulate
for m = 1:NSampPerFrame:NTSample
    sig_idx = m:m+NSampPerFrame-1;
    x1 = jackhammerFileReader();
    x2 = bulldozerFileReader();
    x3 = vibratorFileReader();
    temp = collector([x1 x2],[incidentAngle_jackhammer incidentAngle_bulldozer]);
    sigArray(sig_idx,:) = temp;
    sound_jackhammer(sig_idx) = x1;
    sound_bulldozer(sig_idx) = x2;
    sound_vibrator(sig_idx) = x3;
end

t = 0:1/fs:(size(sigArray, 1)/fs)-1/fs;

%% Locate Sound Source Direction
%Beamform the array.
% mic_loc = getElementPosition(array)';
% 
% lsb = [-80 -80 0];
% usb = [80 80 10];
% 
% NTSample = size(sigArray,1);
% NSamplePerFrame = 10*fs;
% j = 1;
% Searchvolume = {};
% for m = 1:NSamplePerFrame:NTSample
%     if m <= NTSample-NSamplePerFrame+1
%         sig_idx = m:m+NSamplePerFrame-1;
%     elseif m > NTSample-NSamplePerFrame+1
%         sig_idx = m:NTSample;
%     end
%     [finalpos,azimuth_sound,elevation_sound,distance_sound, finalsrp,finalfe,SearchVolume] = srplems(sigArray(sig_idx,:), mic_loc, fs, lsb, usb, c)
%     
%     figure(1);
%     v1=[0,0,0];
%     v2=finalpos;
%     v=[v2;v1];
%     plot3(v(:,1),v(:,2),v(:,3),'b')
%     
%     lsb = [finalpos(1)-5, finalpos(2)-5, finalpos(3)-5];
%     usb = [finalpos(1)+5, finalpos(2)+5, finalpos(3)+5];
%     position(j,:) = finalpos;
%     azimuth(j,:) = azimuth_sound(end);
%     elevation(j,:) = elevation_sound(end);
%     distance(j, :) = distance_sound(end);
%     Searchvolume{j} = SearchVolume;
%     j = j + 1;
% %     scatter3(finalpos(1),finalpos(2),finalpos(3),'filled','MarkerFaceColor','r','MarkerEdgeColor','r');
%     hold on
% end

%% Plot data
% figure(3);
% subplot(6,1,1);
% plot(t,sigArray(:,1))
% title('Signal Received at Channel 1'); ylim([0 duration]);
% axis tight
% 
% subplot(6,1,2);
% plot(t,sigArray(:,2))
% title('Signal Received at Channel 2'); ylim([0 duration]);
% axis tight
% 
% subplot(6,1,3);
% plot(t,sigArray(:,3))
% title('Signal Received at Channel 3'); ylim([0 duration]);
% axis tight
% 
% subplot(6,1,4);
% plot(t,sigArray(:,4))
% title('Signal Received at Channel 4'); ylim([0 duration]);
% axis tight
% 
% subplot(6,1,5);
% plot(t,sigArray(:,5))
% title('Signal Received at Channel 5'); ylim([0 duration]);
% axis tight
% 
% subplot(6,1,6);
% plot(t,sigArray(:,6))
% title('Signal Received at Channel 6'); ylim([0 duration]);
% axis tight

%% 1. Time-delay beamformer
% signalsource = dsp.SignalSource('Signal',sigArray,'SamplesPerFrame',NSampPerFrame);
% TDbeamformer = phased.TimeDelayBeamformer('SensorArray',array,'SampleRate',fs,'Direction',desiredAngle,'PropagationSpeed',c,'WeightsOutputPort',true);
% y_TD = zeros(NTSample,1);
% for m = 1:NSampPerFrame:NTSample
%     [tempTD, w] = TDbeamformer(signalsource());
%     y_TD(m:m+NSampPerFrame-1,:) = tempTD;
% end
% 
% figure(4);
% subplot(2,1,1);
% plot(t,sound_bulldozer)
% xlabel('Time')
% ylabel('Amplitude')
% legend('Original')
% axis tight 
% 
% subplot(2,1,2);
% plot(t,y_TD)
% xlabel('Time')
% ylabel('Amplitude')
% legend('Time-delay beamformer')
% axis tight
% 
% %Calculate the array gain
% agTDB = pow2db(mean((sound_jackhammer+sound_vibrator).^2)/mean((y_TD - sound_bulldozer).^2))

%% 2. Sub-band Phase Shift beamformer
signalsource = dsp.SignalSource('Signal',sigArray,'SamplesPerFrame',NSampPerFrame);
jjj = 1;
subbands = 10:10:500;
% subbands = [25, 50, 100, 150, 200, 250, 300];
% subbands = 295;
for kkk = 1:size(subbands, 2)
    SubbandPhaseShiftbeamformer = phased.SubbandPhaseShiftBeamformer('SensorArray',array, ...
        'Direction',desiredAngle,'OperatingFrequency',oFeq, ...
        'PropagationSpeed',c,'SampleRate',fs,'SubbandsOutputPort',true, ...
        'WeightsOutputPort',true, ...
        'NumSubbands', subbands(kkk));
    reset(signalsource);
    y_SBPS = zeros(NTSample,1);
    for m = 1:NSampPerFrame:NTSample
        [tempSBPS,w,subbandfreq] = SubbandPhaseShiftbeamformer(signalsource());
        y_SBPS(m:m+NSampPerFrame-1,:) = tempSBPS;
    end

%     figure(5);
%     subplot(2,1,1);
%     plot(t,sound_bulldozer)
%     xlabel('Time')
%     ylabel('Amplitude')
%     legend('Original')
%     axis tight
% 
%     subplot(2,1,2);
%     plot(t,real(y_SBPS))
%     xlabel('Time')
%     ylabel('Amplitude')
%     legend('Sub-band phase Shift beamformer')
%     axis tight

    % Calculate the array gain
   agSBPS(jjj, 1) = pow2db(mean((sound_jackhammer+sound_vibrator).^2)/mean((real(y_SBPS) - sound_bulldozer).^2));
%     agSBPS = pow2db(mean((sound_jackhammer+sound_vibrator).^2)/mean((real(y_SBPS) - sound_bulldozer).^2))
   jjj = jjj + 1;
end
plot(subbands, agSBPS ,'r')
title('Array gain vs number of sub-bands','FontSize',26)
xlabel('Number of sub-bands','FontSize',24) 
ylabel('Array gain (dB)','FontSize',24)
set(gca,'fontsize',15)
%% 3. Frost Beamformer
jjj = 1;
filter_lengthFrost = 5:5:100;
% filter_lengthFrost = [5, 10, 20, 40, 60, 80, 100];
for kkk = 1:size(filter_lengthFrost, 2)
    signalsource = dsp.SignalSource('Signal',sigArray,'SamplesPerFrame',NSampPerFrame);
    Frostbeamformer = phased.FrostBeamformer('SensorArray',array,'SampleRate',fs,'PropagationSpeed',c,'FilterLength',filter_lengthFrost(kkk),'DirectionSource','Input port','WeightsOutputPort',true);
    reset(signalsource);
    FrostBeamformer.DiagonalLoadingFactor = 1e-4; % ?????????????????????????????????????
    y_Frost = zeros(NTSample,1);
    for m = 1:NSampPerFrame:NTSample
        [temp, w] = Frostbeamformer(signalsource(),desiredAngle);
        y_Frost(m:m+NSampPerFrame-1,:) = temp;
    end
    
%     figure(6);
%     subplot(2,1,1);
%     plot(t,sound_bulldozer)
%     xlabel('Time')
%     ylabel('Amplitude')
%     legend('Original')
%     axis tight
%     
%     subplot(2,1,2);
%     plot(t,y_Frost)
%     xlabel('Time')
%     ylabel('Amplitude')
%     legend('Frost Beamformer')
%     axis tight
%     
    agFrostB(jjj, 1) = pow2db(mean((sound_jackhammer+sound_vibrator).^2)/mean((y_Frost - sound_bulldozer).^2));
    jjj = jjj + 1;
end
plot(filter_lengthFrost, agFrostB ,'r')
title('Array gain vs filter length','FontSize',26)
xlabel('Filter length','FontSize',24) 
ylabel('Array gain (dB)','FontSize',24)
set(gca,'fontsize',15)
%% 4. Generalized Sidelobe Canceler (GSC) Beamformer
% jjj = 1;
% filter_lengthGSC = [20, 40, 60, 80, 100]; % 55
% for kkk = 1:size(filter_lengthGSC, 2)
%     signalsource = dsp.SignalSource('Signal',sigArray, 'SamplesPerFrame',NSampPerFrame);
%     LMSStepSize = 0.1;
%     GSCbeamformer = phased.GSCBeamformer('SensorArray',array,'PropagationSpeed',...
%         c,'SampleRate',fs,'Direction',desiredAngle,...
%         'FilterLength',filter_lengthGSC(kkk), 'LMSStepSize', LMSStepSize);
%     reset(signalsource);
%     
%     y_GSC = zeros(NTSample,1);
%     for m = 1:NSampPerFrame:NTSample
%         temp = GSCbeamformer(signalsource());
%         y_GSC(m:m+NSampPerFrame-1,:) = temp;
%     end
%     
%     % figure(7);
%     % subplot(2,1,1);
%     % plot(t,sound_bulldozer)
%     % xlabel('Time')
%     % ylabel('Amplitude')
%     % legend('Original')
%     %
%     % subplot(2,1,2);
%     % plot(t,y_GSC)
%     % xlabel('Time')
%     % ylabel('Amplitude')
%     % legend('GSC Beamformer')
%     
%     % Calculate the array gain
%     agGSCB(jjj, 1) = pow2db(mean((sound_jackhammer+sound_vibrator).^2)/mean((y_GSC - sound_bulldozer).^2));
%     jjj = jjj + 1;
% end
% plot(filter_lengthGSC, agGSCB)
%% 5. Time Delay Linear Constraint Minimum Variance (LCMV) Beamformer
% jjj = 1;
% filter_lengthLCMV = [20, 40, 60, 80, 100]; % 30
% for kkk = 1:size(filter_lengthLCMV, 2)
%     constraintMatrix = kron(eye(filter_lengthLCMV(kkk)),ones(Nele,1));
%     desiredResponseVector = eye(filter_lengthLCMV(kkk),1);
%     signalsource = dsp.SignalSource('Signal',sigArray, 'SamplesPerFrame',NSampPerFrame);
%     TimeDelayLCMVbeamformer = phased.TimeDelayLCMVBeamformer('SensorArray',array,...
%         'PropagationSpeed',c,'SampleRate',fs,'FilterLength',filter_lengthLCMV(kkk),...
%         'Direction',desiredAngle,'Constraint',constraintMatrix,...
%         'DesiredResponse',desiredResponseVector);
%     reset(signalsource);
%     
%     y_TDLCMV = zeros(NTSample,1);
%     for m = 1:NSampPerFrame:NTSample
%         temp = TimeDelayLCMVbeamformer(signalsource());
%         y_TDLCMV(m:m+NSampPerFrame-1,:) = temp;
%     end
%     
%     % figure(8);
%     % subplot(2,1,1);
%     % plot(t,sound_bulldozer)
%     % xlabel('Time')
%     % ylabel('Amplitude')
%     % legend('Original')
%     %
%     % subplot(2,1,2);
%     % plot(t,y_TDLCMV)
%     % xlabel('Time')
%     % ylabel('Amplitude')
%     % legend('TDLCMV Beamformer')
%     
%     % Calculate the array gain
%     agTDLCMV(jjj, 1) = pow2db(mean((sound_jackhammer+sound_vibrator).^2)/mean((y_TDLCMV - sound_bulldozer).^2));
%     jjj = jjj + 1;
% end
% plot(filter_lengthLCMV, agTDLCMV)
%% 6. Subband MVDR Beamformer
jjj = 1;
NumSubbands = [5 ,10, 15, 20, 25, 30, 35]; % 20
for kkk = 1:size(NumSubbands, 2)
    signalsource = dsp.SignalSource('Signal',sigArray, 'SamplesPerFrame',NSampPerFrame);
    SubbandMVDRbeamformer = phased.SubbandMVDRBeamformer('SensorArray',array,...
        'Direction',desiredAngle,'OperatingFrequency',oFeq,...
        'PropagationSpeed',c,'SampleRate',fs, ...
        'SubbandsOutputPort',true,'WeightsOutputPort',true, 'NumSubbands', NumSubbands(kkk));
    reset(signalsource);
    
    y_SBMVDR = zeros(NTSample,1);
    for m = 1:NSampPerFrame:NTSample
        [temp,w,subbandfreq] = SubbandMVDRbeamformer(signalsource());
        y_SBMVDR(m:m+NSampPerFrame-1,:) = temp;
    end
    
%     figure(9);
%     subplot(2,1,1);
%     plot(t,sound_bulldozer)
%     xlabel('Time')
%     ylabel('Amplitude')
%     legend('Original')
%     axis tight
%     
%     subplot(2,1,2);
%     plot(t,real(y_SBMVDR))
%     xlabel('Time')
%     ylabel('Amplitude')
%     legend('TDLCMV Beamformer')
%     axis tight
%     
    % Calculate the array gain
    agSBMVDR(jjj, 1) = pow2db(mean((sound_jackhammer+sound_vibrator).^2)/mean((real(y_SBMVDR) - sound_bulldozer).^2));
    jjj = jjj + 1;
end
plot(NumSubbands, agSBMVDR)
