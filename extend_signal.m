function sigExtended = extend_signal(sigMatrix, extFactor)
% extendsignal - extend the emg signals with given extension factor
% input:
%       sigMatrix: the high-density electromyography in channels x time
%       extFactor: the extension factor
% output:
%       sigExtended: extended signal
% 
% Author: Yue Wen at University of Central Florida, Nov. 20, 2023.
%

% set default extension factor if not provided
if nargin == 1
    % extFactor = 16;
    extFactor = round(1024/min(size(sigMatrix)));
end

% get signal channels and length
[nFeatures, nSamples] = size(sigMatrix);
if nSamples < nFeatures
    sigMatrix = sigMatrix';
    [nFeatures, nSamples] = size(sigMatrix);
end

% aligned at the beginning
signal = [zeros(nFeatures, extFactor-1), sigMatrix];
% aglined at the center
% signal = [zeros(nfeatures, extfactor/2), SIG_m, zeros(nfeatures, extfactor/2)];
signalExtended = zeros(extFactor, nFeatures, nSamples);
for iSample = 1:nSamples
    signalExtended(:, :, iSample) = signal(:, iSample:iSample+extFactor-1)';
end
sigExtended = reshape(signalExtended, [nFeatures*extFactor, nSamples]);
end