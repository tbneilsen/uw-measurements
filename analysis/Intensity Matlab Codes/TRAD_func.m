givefunction TRAD = TRAD_func(fss,Xss,probe_config,rho,c)
%{
This function uses the traditional method to calculate estimated intensity
vectors from multi-microphone intensity probes given single-sided ffts

TRAD - short for traditional method

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


OUTPUT: a structure, TRAD, estimating acoustic intensity. It contains the 
following variables:

I - Active time averaged intensity

Q - Reactive time averaged intensity

Ep - Potential energy density (same for TRAD and PAGE)

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

omega = 2*pi.*fss;

size_fft = size(Xss);
numblocks = size_fft(2);
numCH = size_fft(1);
numsamp = 2*size_fft(3);



%%% The following calculates X and the pseudo inverse of X, used in
%%% calculating a least-squares estimate of grad_p.
%%% See write-up for more details on this and the following processing
%%% techniques

%Calculates X and the pseudo inverse of X
numdif = (numCH^2-numCH)/2; %Number of unique microphone pairs
index = 1;
X = zeros(numdif,3); %Array of distances between microphone pairs.
for ii = 1:numCH
    for jj = ii+1:numCH
        X(index,:) = probe_config(jj,:)-probe_config(ii,:);
        index = index+1;
    end
end
pinvX = pinv(X); %pseudo inverse of X

%%%Initializing of variables
I_blocks = zeros(numblocks,3,numsamp/2);
Q_blocks = zeros(numblocks,3,numsamp/2);
psq_blocks = zeros(numblocks,numsamp/2);
usq_blocks = zeros(numblocks,3,numsamp/2);
pdiff = zeros(numdif,numsamp/2);

%%% The following calculates the TRAD intensities of each block. This is
%%% done by calculating grad_p and p0 for each block.
for m = 1:numblocks
    p = squeeze(Xss(:,m,:));

    %grad_p
    index = 1;
    for ii = 1:numCH
        for jj = ii+1:numCH
            pdiff(index,:) = p(jj,:) - p(ii,:);%calculates the differences in complex pressures - for the TRAD method
            index = index+1;
        end
    end
    
    grad_p = pinvX*pdiff;
    
    %p0
    centerCH = find((probe_config(:,1)==0)&(probe_config(:,2)==0)&(probe_config(:,3)==0),1);
    if length(centerCH)==1;
        %If there is a microphone at the center of the probe
        p0 = p(centerCH,:); %complex p0
        P0 = abs(p(centerCH,:));%abs p0
    else
        %If there is NOT a microphone at the center of the probe
        p0 = mean(p,1); 
        P0 = mean(abs(p),1);
    end
    
    %The function repmat repeats omega and p0 across x, y and z so they can
    %be multiplied by grad_p element-wise in the intensity calculation.
    omega_vec = repmat(omega,3,1);
    p0_vec = repmat(p0,3,1);
    
    %Intensity for the current block
    Ic = conj(1j./(omega_vec*rho).*grad_p).*p0_vec; %Complex intensity
    %Ic = 1j./(omega_vec*rho).*conj(grad_p).*p0_vec; %Complex intensity
    u = 1j./(omega_vec*rho).*grad_p; %Particle velocity
    
    I_blocks(m,:,:) = real(Ic); %Real part of intensity
    Q_blocks(m,:,:) = imag(Ic); %Imaginary part of intensity
    
    usq_blocks(m,:,:) = abs(u).^2; %Define u^2 and p^2 for energy quantities
    psq_blocks(m,:) = P0.^2;%abs(p0).^2;%Accurate at higher frequencies if P0 instead of p0

end

%Averages the intensities across blocks
I = squeeze(mean(I_blocks,1));
Q = squeeze(mean(Q_blocks,1));

%Calculates energy densities
psq = mean(psq_blocks,1);
usq = squeeze(mean(usq_blocks,1));
Ep = psq./(2*rho*c^2); %Potential Energy Density
Ekvec = rho/2*usq; 
Ek= sum(Ekvec,1); %Kinetic Energy Density
E = Ep+Ek; %Total Energy Density
EL = Ek - Ep; %Lagrangian Energy Density 

%Calculates magnitude of intensity vector
I_mag = squeeze(sqrt(real(I(1,:)).^2 +real(I(2,:)).^2+real(I(3,:)).^2));
Q_mag = squeeze(sqrt(real(Q(1,:)).^2 +real(Q(2,:)).^2+real(Q(3,:)).^2));

%Finds the angle (in radians) of the intensity vector (in x-y plane)
Ix = I(1,:);
Iy = I(2,:);
I_dir = atan2(Iy,Ix);
Qx = Q(1,:);
Qy = Q(2,:);
Q_dir = atan2(Qy,Qx);

%Creates structure with all the variables
TRAD.I = I;
TRAD.Q = Q;
TRAD.Ep = Ep;
TRAD.Ek = Ek;
TRAD.E = E;
TRAD.EL = EL;
TRAD.I_mag = I_mag; 
TRAD.I_dir = I_dir;
TRAD.Q_mag = Q_mag;
TRAD.Q_dir = Q_dir;
TRAD.p0 = p0;
TRAD.grad_p = grad_p;
TRAD.u = u;
TRAD.z=repmat(p0,[3,1])./u;