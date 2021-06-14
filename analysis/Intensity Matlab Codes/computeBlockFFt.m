function [Xss,numblocks] = computeBlockFFt(data,ns,w,W)
% INPUT:
% data is the time domain data. It should be a single vector.
% ns is the blocksize
% w is the window
% W is the weight (W = mean(W.^2))
% OUTPUT:
% Xss is the single sided scaled data
% numblocks is the number of blocks that the data was broken up into

% Make sure that the is in the right format.
if size(data,1) == 1 && size(data,2) > 1
    data = data';
end

%%% ZERO MEAN
data = data-mean(data);
%%% LENGTH OF DATA
N = length(data);
%%% NUMBER OF BLOCKS
numblocks = floor(2*N/ns-1);
%%% CHUNK DATA
blockData = [reshape(data(1:((numblocks)*ns/2)),[ns/2,(numblocks)])',...
    reshape(data((1:((numblocks)*ns/2))+ns/2),[ns/2,numblocks])'];
blockData = blockData.*repmat(w,[numblocks,1]);
%%% FFT
X = fft(blockData,ns,2);

%%% Scale data so that conj(Xss).*Xss = Gxx
Xss = X(:,1:floor(ns/2))/ns*sqrt(2/W);

end