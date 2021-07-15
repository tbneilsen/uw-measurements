lims = [100 1000];
t = 0:.005:100;
fs = 48000;
data = sin(t);
[fc, OTOspec] = TDOTOspec(data, fs, lims)