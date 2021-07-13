function PAGE = PAGE_func(fss,Xss,probe_config,rho,c,opt,unwrap_num)
%{
This function uses the PAGE method to calculate estimated intensity vectors
from multi-microphone intensity probes given single-sided ffts

PAGE - short for Phase and Amplitude Gradient Estimates

Equation numbers referenced in the code are from "Phase and amplitude
gradient method for the estimation of acoustic vector quantities" (Thomas
et. al., JASA, Jun 2015)

INPUTS:

fss - single-sided frequency array corresponding to Xss

Xss - matrix of blocked, single-sided ffts, organized as:
    [CH_number,block_number,sample]
For example, the fft of the 2nd microphone's 15th block would be:
    Xss(2,15,:)
The ffts should already be scaled correctly so that when two ffts are 
multiplied, the result is a cross-spectrum. The order of the channels 
should correspond to the configuration variable: probe_config

probe_config - An nx3 matrix, where n is the number of channels per probe. 
It contains the location of each microphone relative to the geometric
center of the probe, in [x,y,z]. If there is a center microphone, it should
have the coordinates [0,0,0].

rho - density of air

c - speed of sound

opt - Option that determines which unwrapping method is used by unwrap_func.
Possible values:
0 - no unwrapping
1 - use MATLAB's built-in unwrapping function (default)
average - determines unwrapping by the average of the last several points
slope - determines unwrapping by the slope of the last several points

unwrap_num - Number of points to consider in the average and slope
unwrapping methods. Defaults to 8.


OUTPUT: a structure, PAGE, estimating acoustic intensity. It contains the 
following variables:

I - Active time averaged intensity

Q - Reactive time averaged intensity

Ep - Potential energy density (same for FD and PAGE)

Ek - Kinetic energy density

E - Total energy density (Ep + Ek)

EL - Lagrangian energy density (Ek - Ep)
%}


if nargin < 4
    rho = 1.21;
end

if nargin < 5
    c = 343;
end

% If a value is not supplied for opt, it defaults to unwrapping
if nargin < 6
    opt = 1;
end

if nargin < 7
    unwrap_num = [];
end

omega = 2*pi.*fss;

size_fft = size(Xss);
numCH = size_fft(1);
numsamp = 2*size_fft(3);


%%% The following quantities are calculated for estimating intensity:
%%%     grad_phi (gradient of the phase)
%%%     P0 (center pressure amplitude)
%%%     grad_P (gradient of the pressure amplitude)
%%% (Each of these quantities are frequency dependent and are represented
%%% by a vector)


%%% grad_phi

%%% First, grad_phi is calculated by a least-squares estimate. This requires
%%% the pseudo inverse of X, where X is an array of distances between
%%% microphone pairs. This is multiplied by phidiff, which contains the
%%% differences in phase between each pair as calculated by their transfer
%%% function.
%%% See write-up for more details on this and the following processing
%%% techniques.

%Calculates X and the pseudo inverse of X
numdif = (numCH^2-numCH)/2; %Number of unique microphone pairs
index = 1;
X = zeros(numdif,3);
for ii = 1:numCH
    for jj = ii+1:numCH
        X(index,:) = probe_config(jj,:)-probe_config(ii,:); %(Eq. 2)
        index = index+1;
    end
end
pinvX = pinv(X); %pseudo inverse of X

%Calculates the transfer function of each microphone pair as <Gxy>/<Gxx>
Gxx = mean(conj(Xss).*Xss,2);
TF = zeros(numdif,numsamp/2);
coh = zeros(numdif,numsamp/2);
count = 1;
for ii = 1:numCH
    for jj = (ii+1):numCH
        Gxy = squeeze(mean(conj(Xss(ii,:,:)).*Xss(jj,:,:),2)).';
        TF(count,:) = Gxy./Gxx(ii,:);
        coh(count,:) = abs(Gxy).^2./Gxx(ii,:)./Gxx(jj,:);
        count = count + 1;
    end
end

%Calculates the spatial Nyquist frequency of the probe
fN = c/max(squeeze(sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2)))/2;

%Calculates phidiff as the phase of the transfer function
%Unwrapping applied at the Nyquist frequency
unwraplim = find(fss>=fN*0.9,1); %Find index to start unwrapping
phidiff = unwrap_func(-angle(TF),[],2,opt,unwrap_num,coh,unwraplim); %(Eq. 16)

%Finds the least-squares estimate of the phase gradient
grad_phi = pinvX*phidiff; %(Eq. 12)


%%% P0

%%% Next, the center pressure is estimated for each block, then averaged
%%% over the blocks. If a center microphone exists, it is used for the
%%% center pressure, otherwise the rms value of all microphones is used.

%Calculate center pressure for each block
centerCH = find((probe_config(:,1)==0)&(probe_config(:,2)==0)&(probe_config(:,3)==0),1);
if length(centerCH)==1;
    %If there is a microphone at the center of the probe
    P0_blocks = abs(Xss(centerCH,:,:)); %Use center pressure amplitude
else
    %If there is NOT a microphone at the center of the probe
    P0_blocks = mean(abs(Xss),1);
end

%Averages center pressure across blocks
P0 = squeeze(sqrt(mean(P0_blocks.^2,2)))';


%%% grad_P

%%% To calculate grad_P, we need the pressure amplitude at each microphone.
%%% The pseudo inverse of X is multiplied by an array of the differences of
%%% these pressure amplitudes.

%Finds the magnitude of the pressures for each block and channel
P_blocks = abs(Xss); %(Eq. 18)

%Averages the pressure amplitudes across blocks for each channel
P = squeeze(sqrt(mean(P_blocks.^2,2)));

%Finds the least-squares estimate of the amplitude gradient
Pdiff = zeros(numdif,numsamp/2);
index = 1;
for ii = 1:numCH
    for jj = ii+1:numCH
        Pdiff(index,:) = P(jj,:) - P(ii,:); %(Eq. 20)
        index = index+1;
    end
end
grad_P = pinvX*Pdiff; %(Eq. 19)


%%% Now that we have the needed physical quantities, we calculate the
%%% intensities and energy densities, and store them to a structure to output.

%The function repmat repeats omega and P0 across x, y and z so they can be
%multiplied by grad_phi and grad_P element-wise.
omega_vec = repmat(omega,3,1);
P0_vec = repmat(P0,3,1);

%Calculates the intensities (no 1/2 because these are rms pressures)
I = 1./(omega_vec*rho).*(P0_vec.^2).*grad_phi; %Active intensity (Eq. 11a)
Q = -1./(omega_vec*rho).*P0_vec.*grad_P; %Reactive intensity (Eq. 11b)
u = 1./(omega_vec*rho).*(P0_vec.*grad_phi+1j*grad_P); %for energy quantities

%Calculates energy densities
Ep = (P0.^2)./(2*rho*c^2); %Potential Energy Density
Ekvec = rho/2*real(u.*conj(u));
Ekvec(:,1) = 0;
Ek = sum(Ekvec,1); %Kinetic Energy Density
E = Ep+Ek; %Total Energy Density
EL = Ek - Ep; %Lagrangian Energy Density

%Calculates magnitude of intensity vector
I_mag = squeeze(sqrt(real(I(1,:)).^2 +real(I(2,:)).^2+real(I(3,:)).^2));
I_mag(1) = 0;
Q_mag = squeeze(sqrt(real(Q(1,:)).^2 +real(Q(2,:)).^2+real(Q(3,:)).^2));
Q_mag(1) = 0;

%Finds the angle (in radians) of the intensity vector (in x-y plane)
Ix = I(1,:);
Iy = I(2,:);
I_dir = atan2(Iy,Ix);
Qx = Q(1,:);
Qy = Q(2,:);
Q_dir = atan2(Qy,Qx);

%Creates structure with all the variables
PAGE.I = I;
PAGE.Q = Q;
PAGE.Ep = Ep;
PAGE.Ek = Ek;
PAGE.E = E;
PAGE.EL = EL;
PAGE.I_mag = I_mag;
PAGE.I_dir = I_dir;
PAGE.Q_mag = Q_mag;
PAGE.Q_dir = Q_dir;
PAGE.fN = fN;
PAGE.TFphase = phidiff;
PAGE.TFphase_wrapped = -angle(TF);
PAGE.P0 = P0;
PAGE.grad_P = grad_P;
PAGE.grad_phi = grad_phi;
PAGE.u = u;
PAGE.z=repmat(P0.^2,[3,1])./(I+1j*Q);