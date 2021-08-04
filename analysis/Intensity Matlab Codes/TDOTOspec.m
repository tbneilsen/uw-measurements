function [ fc,OTOspec ] = TDOTOspec(x,fs,flims)
% OTO filter calculations for time-domain filtering
% Kent L. Gee, kentgee@byu.edu
% Modified 24 Oct. 2013

% Define preferred band center frequencies
fc = [                                       1   1.25    1.6    2.0...
    2.5  3.15    4   5   6.3   8     10   12.5    16    20 ...
    25   31.5   40   50   63   80   100   125     160   200 ...
    250  315    400  500  630  800  1000  1250    1600  2000 ...
    2500 3150   4000 5000 6300 8000 10000 12500   16000 20000 ...
    25e3 31.5e3 40e3 50e3 63e3 80e3 100e3];


% Truncate fc array
ind = find(fc>=flims(1) & fc<=flims(2));
fc=fc(ind);

for n=1:length(fc)
    if fc(n)<fs/20  %factor of 20 suggested as threshold for decimation
        decirat(n)=ceil(fs/fc(n)/20); %decimation factor
    else
        decirat(n)=1;
    end
    fsdec(n)=fs/decirat(n);        %decimated sampling rate
    
    % Design Butterworth 2Nth-order one-third-octave filter, for order N=3
    % Note: BUTTER is based on a bilinear transformation, as suggested in
    % ANSI standard.  From oct3dsgn function by Christophe Couvreur, Faculte Polytechnique de Mons (Belgium)
    N=3;
    f1 = fc(n)/(2^(1/6));
    f2 = fc(n)*(2^(1/6));
    Qr = fc(n)/(f2-f1);
    Qd = (pi/2/N)/(sin(pi/2/N))*Qr;
    alpha = (1 + sqrt(1+4*Qd^2))/2/Qd;
    W1 = fc(n)/(fsdec(n)/2)/alpha;
    W2 = fc(n)/(fsdec(n)/2)*alpha;
    [B(n,:),A(n,:)] = butter(N,[W1,W2]);
    
end

% Zero-mean, window, and scaling of data by equivalent noise bandwidth

x=x-mean(x);
win=hanning(length(x))';
x=x.*win/sqrt(mean(win.^2));


for n=1:length(fc)
    if decirat(n)>1
        xdec=resample(x,1,decirat(n));
    else
        xdec=x;
    end
   % explored prepending data, but with windowing, doesn't matter. 
   % Nd=length(xdec);
   % xfilt=filter(B(n,:),A(n,:),[xdec,-xdec(end:-1:1)]);                 
   % OTOspec(n)=mean(xfilt(floor(Nd/2)+1:floor(Nd/2)+1+Nd).^2); %calulate mean-square value for non-prepended block
    
    xfilt=filter(B(n,:),A(n,:),[xdec]);               %filter waveform
    OTOspec(n)=mean(xfilt.^2); %calulate mean-square
end

end

