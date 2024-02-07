function [signalMatrix, signalExtended, signalWhitended, whitenParameters] = extend_and_whiten(signalCell, discardChannelsVec)
%% extend and whiten the EMG signals
% input: 
%       signalCell: mxn cell corresponding to array electrodes
%       discardChannelsVec: mxn matrix indicting the discarded channels
% output:
%       signalMatrix: filtered SIG in matrix format
%       signalExtended: filtered and extended SIG
%       signalWhitended: filtered, extended, and whiened SIG
% 
% Author: Yue Wen at University of Central Florida, Nov. 20, 2023.
%

% required parameters
differentialMode = 0;
nbextchan = 1000;
notchFrequency = 50;
sampleRate = 2048;

% SIG (13*5 cell) -> EMGs (65*n matrix)
iChannel = 1;
for r=1:size(signalCell,1)
    for c=1:size(signalCell,2)
        if ~isempty(signalCell{r,c}) && discardChannelsVec(r,c) == 0
            signalMatrix(iChannel,:) = signalCell{r,c};
            iChannel = iChannel+1;
        end
    end
end

% Removing the mean from each channel
meanValue = mean(signalMatrix, 2);
signalMatrix = signalMatrix - meanValue;

%% Step 1: Filter signal
% Filter signal with a Bandpass filter: preserve signals in range of 20 Hz and 500 Hz
[b,a] = butter(2, [20 500]/(sampleRate/2)); % Bandpass Butterworth 2nd order, Cutoff 10-500Hz
for j=1:size(signalMatrix,1)
    signalMatrix(j,:) = filtfilt(b, a, signalMatrix(j,:)); % Zero-phase digital filtering
end

% Remove line interference with a Notch filter
wo = notchFrequency/(sampleRate/2);
bw = wo/35;         % Q factor for iirnotch filter
[b, a] = iirnotch(wo, bw);
signalMatrix = filter(b, a, signalMatrix, [], 2);

% Differentiation mode
if differentialMode == 1
    signalMatrix = diff(signalMatrix, 1, 2);
end

%% Step 2: Extend signal
%  Extend signal with extension factor to reach 1000+ channels
extendFactor = round(nbextchan/size(signalMatrix,1));
signalExtended = extend_signal(signalMatrix, extendFactor);

%% Step 3: Whiten signal
% Whitening with ZCA
[signalWhitended, whiteningMatrix, dewhiteningMatrix] = whiten_signal(signalExtended);

% Remove the edges
signalMatrix = signalMatrix(:,1:end);
signalExtended = signalExtended(:,1:end);
signalWhitended = signalWhitended(:,1:end);
whitenParameters.whiteningMatrix = whiteningMatrix;
whitenParameters.dewhiteningMatrix = dewhiteningMatrix;
end