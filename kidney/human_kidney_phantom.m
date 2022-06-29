% Creates a phantom for a liver from a MR scan of the liver.
% The size of the image is 100 x 100 (width and depth) 
% and the thickness is 15 mm. 
% The phantom starts 2 mm from the transducer surface.
% 
% Ver. 1.1, March 29, 2000, Jørgen Arendt Jensen

function [positions, amp] = human_kidney_phantom (N)

% Load the bitmap image

[liv_kid, MAP]=bmpread('kidney_cut.bmp');

% Define image coordinates

liv_kid=liv_kid';
[Nl, Ml]=size(liv_kid);

x_size = 100/1000 ;    %  Size in x-direction [m]
dx=x_size/Nl;          %  Sampling interval in x direction [m]
z_size = 100/1000 ;    %  Size in z-direction [m]
dz=z_size/Ml;          %  Sampling interval in z direction [m]
y_size = 15/1000;      %  Size in y-direction [m]
theta = 35/180*pi;     %  Rotation of the scatterers [rad]
theta = 0;
z_start = 2/1000;


% Calculate position data

x0 = rand(N,1);
x = (x0-0.5)* x_size;
z0 = rand(N,1);
z = z0*z_size+z_start;
y0 = rand(N,1);
y = (y0-0.5)* y_size; 

%  Find the index for the amplitude value

xindex = round((x + 0.5*x_size)/dx + 1);
zindex = round((z - z_start)/dz + 1);
inside = (0 < xindex)  & (xindex <= Nl) & (0 < zindex)  & (zindex <= Ml);
index = (xindex + (zindex-1)*Nl).*inside + 1*(1-inside);

% Amplitudes with different variance must be generated according to the the 
% input map.
% The amplitude of the liver-kidney image is used to scale the variance

amp=exp(liv_kid(index)/100);
amp=amp-min(min(amp));
amp=1e6*amp/max(max(amp));
amp=amp.*randn(N,1).*inside;

%  Generate the rotated and offset block of sample

xnew=x*cos(theta)+z*sin(theta);
znew=z*cos(theta)-x*sin(theta);
znew=znew-min(min(znew)) + z_start;

positions=[(xnew-40/1000) y znew];
positions=[xnew y znew];


