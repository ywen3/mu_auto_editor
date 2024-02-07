function [signalWhitened, whiteningMatrix, dewhiteningMatrix] = whiten_signal(signalExtended, varianceThreshold)
% whiten_signal: (ZCA) whiten the high-density electromyography
% input:
%       signalExtended: extended high-density electromyography
% output:
%       signalWhitened: whitened signal
%       whiteningMatrix: 
%       dewhiteningMatrix:
% 
% Author: Yue Wen at University of Central Florida, Nov. 20, 2023.
% 

% pca components
if nargin == 1
    varianceThreshold = 1;
end
[E, D]= pcamat(signalExtended);
% [E, ~, D] = pca(signalExtended);

% Get a subset of components
if varianceThreshold < 1
    eigenvalues = diag(D);
    [dvalue, dindex] = sort(eigenvalues, 'descend');
    dvalue = dvalue./sum(dvalue);
    d = 0;
    k = 1;
    while d<varianceThreshold
        d = d + dvalue(k);
        k = k + 1;
    end
    subindex = dindex(1:k);
    eigenVectors = E(:, subindex);
    eigenValueMatrix = diag(eigenvalues(subindex));
else
    eigenVectors = E;
    eigenValueMatrix = D;
end

% caclulate whitening and de-whitenning matrix
sqrtEigenValueMatrix = sqrt(eigenValueMatrix);
whiteningMatrix = eigenVectors*inv(sqrtEigenValueMatrix)*eigenVectors';
dewhiteningMatrix = eigenVectors*sqrtEigenValueMatrix*eigenVectors';
signalWhitened = whiteningMatrix*signalExtended;

end