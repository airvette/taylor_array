function [u, theta, arrayFactor] = TaylorSynth(f0, theta0, theta3db, slldb)
% Models the relative power levels of a Taylor Illuminated Phased Array
%  based on frequency, scanning direction, the beamwidth and the Sidelobe
%  Level.  References for this function are Bro08 (Chapter 7), Mai05
%  (Chapter 3), and Sko01 (Chapter 9.5-9.16)
%================ Inputs =========================================
% f0 = the operating frequency of the array
% theta0 = the scan angle of the array, in degrees.  This is the angle off
%  of the boresight in which the main beam is pointing
% theta3db = the angle off of boresight where the mainbeam is 3dB less than
%  the maximum value.  Also known as the beamwidth.  This value is entered
%  in degrees
% slldb = the sidelobe level, which indicates the height of the first side
%  lobe relative to the main beam.  This value is in decibles and is <0.
%================ Outputs ========================================
% u = the angle range the values are computed over.  Proportional to sin(theta)
% theta = the angle range the valuse are computed over
% arrayFactor = the Taylor Illuminated power values
%////////////////////////////////////////////////////////////////
% Jeffrey Guido, UCCS Masters Thesis, FA 2013
 
theta = -pi/2:pi/10000:pi/2; % initialitize angle range
c = 3e8; % speed of light, in m/s
lambda0 = c/f0; % calculation of the wavelength, in m
 
% Determine the Taylor Parameter, A
A = (1/pi)*acosh(10^(-slldb/20)); %Bro08, pg53, eq(7.5)
 
% Determine the Taylor Order, nBar
nBar = round((slldb/22.8)^2 - (slldb/36.3) + .759); % Bro08, pg63, eq(7.32)
 
% Limit the size of nBar
if nBar > 1000
    nBar = 1000;
end
 
% Determine the Taylor Parameter, sigmaSq
sigmaSq = (nBar^2)/(A^2 + (nBar-.5)^2); %Bro08, pg53, eq(7.3)
 
% Determine the Taylor Coefficients, fVector
fVector = ones(nBar, 1); % There are nBar - 1 Taylor Coefficients plus the 
    % zero-th coefficient which is equal to 1
for k = 1:length(fVector)-1 % calc the Taylor Coeff's F1 through F(nBar-1)
    numProd = 1; % starting value for the numerator
    denProd = 2; % starting value for the denominator 
    for n = 1:nBar-1 % iterate through the numerator and denominator products
        % Calc numerator array product
        numProd = numProd * (1-(k^2/(sigmaSq*(A^2+(n-.5)^2)))); 
        if n ~= k
            % Calc denominator array product
            denProd = denProd * (1-(k^2/n^2)); 
        end
    end   
    % Perform final Taylor Coeff calculation for this zero
    fVector(k+1) = (-1)^(k+1)*numProd/denProd;
end    
 
% WRITE CODE FOR FINDING THE F VECTOR
%fVector = [1; .387482; -.00956429; .0046963; -.00133399];
 
u3db = HalfPowerPt(fVector); % Determine the value of u where the beam is 
    % at its half power
 
a = (lambda0*u3db)/(2*sind(theta3db)); % compute half-width of array
 
u = (2*a/lambda0).*(sin(theta)-sind(theta0));
 
 
% Prep vectors for a vectorized calculation
% Prep F Vector for calculation
reflectedFVector = [zeros(length(fVector)-1,1);fVector]; % allocate space
for i = 2:length(fVector)
    reflectedFVector(length(fVector)-i+1) = fVector(i); % write values to 
        % the reflected vector
end
 
mVec = -(length(fVector)-1):(length(fVector)-1); % create a vector of m values
    % for the below while loop, ref Bro08, pg 71.
arrayFactor = zeros(length(u),1); % allocate space for output
u_tmp = zeros(length(u),length(mVec)); % allocate space for operational array
% Conduct the array factor calculation for each value of u
for j = 1:length(u)
    u_tmp(j,:) = u(j); % populate array rows
    u_tmp(j,:) = u_tmp(j,:) + mVec; % add in the m values
    FuComponents = sinc(u_tmp(j,:)).*reflectedFVector';
    arrayFactor(j) = sum(FuComponents); % sum the calculations and return the 
        % array factor
end
 
function halfPwr = HalfPowerPt(fVector)
% Determines the point of half power in terms of u, the value of sin(theta)
%  off boresight.  
 
u3db = 0; % dummy starting value for the half power iteration
u3db_next = .55; % starting variable for the next iteration
 
% Prep F Vector for calculation
reflectedFVector = [zeros(length(fVector)-1,1);fVector]; % allocate space
for i = 2:length(fVector)
    reflectedFVector(length(fVector)-i+1) = fVector(i); % write values to 
        % the reflected vector
end
 
mVec = -(length(fVector)-1):(length(fVector)-1); % create a vector of m values 
    % for the below while loop, ref Bro08, pg 71.
 
while abs(u3db_next - u3db) > .001 % iterate until u3db converges
    u3db = u3db_next; % pass the value found last last iteration to the current value
    % Find the Array Factor Value at the Projected Half Power Point
    fU3db = sum(reflectedFVector.*((sin(pi.*(u3db+mVec'))./(pi.*(u3db+mVec')))));
    % Find the Derivative of the Array Factor
    % Derivative Numerator
    delFUnum = (pi.*(u3db+mVec').*cos(pi.*(u3db+mVec'))-sin(pi.*(u3db+mVec'))); 
    % Derivative Denominator
    delFUden = pi.*(u3db+mVec').^2;
    % Complete the Derivative Calculation
    delFU3db = sum(reflectedFVector.*(delFUnum./delFUden));
    % Find the next u3db value
    u3db_next = u3db - (fU3db-(1/sqrt(2)))/delFU3db;
end
 
halfPwr = u3db_next; % return the half power point
