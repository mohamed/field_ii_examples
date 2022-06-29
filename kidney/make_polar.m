%  Make the image interpolation for the polar scan
%
%  Version 1.0, 17/2-99 JAJ
%  Version 1.1, 16/8-2007 JAJ
%    Small changes in compression

%  Set initial parameters

D=10;       %  Sampling frequency decimation factor
fs=100e6/D; %  Sampling frequency  [Hz]
c=1540;     %  Speed of sound [m]
no_lines=128;                  %  Number of lines in image
image_width=90/180*pi;         %  Size of image sector [rad]
dtheta=image_width/no_lines;   %  Increment for image

%  Read the data and adjust it in time 

min_sample=0;
for i=1:no_lines

  %  Load the result

  cmd=['load rf_data/rf_ln',num2str(i),'.mat']
  eval(cmd)
  
  %  Find the envelope
  
  rf_env=abs(hilbert([zeros(round(tstart*fs-min_sample),1); rf_data]));
  env(1:max(size(rf_env)),i)=rf_env;
  end

%  Do logarithmic compression to 40 dB

dB_Range=50;
env=env-min(min(env));
log_env=20*log10(env(1:D:max(size(env)),:)/max(max(env)));
log_env=255/dB_Range*(log_env+dB_Range);

%  Get the data into the proper format

start_depth=0.02;   % Depth for start of image in meters
image_size=0.105;   % Size of image in meters
skipped_samples=0;
samples=max(size(log_env));

start_of_data=(skipped_samples+1)/fs*c/2;           % Deepth for start of data in meters
end_of_data=(skipped_samples+samples+1)/fs*c/2;     % Depth for end of data in meters
delta_r=c/2*1/fs;                                   % Sampling interval for data in meters

theta_start= -no_lines/2*dtheta;     % Angle for first line in image

Nz=512;                         % Size of image in pixels
Nx=512;                         % Size of image in pixels
scaling=1;                      % Scaling factor form envelope to image

[N,M]=size(log_env);
D=floor(N/1024);
env_disp=uint8(255*log_env(1:D:N,:)/max(max(log_env)));

%  Make the tables for the interpolation

make_tables (start_depth, image_size,          ...
             start_of_data, delta_r*D, round(samples/D),  ...
	     -theta_start, -dtheta, no_lines,       ...
	     scaling, Nz, Nx);

%  Perform the interpolation

img_data=make_interpolation (env_disp);

%  Display the image 

image(((1:Nz)/Nz-0.5)*image_size*1000, ((1:Nx)/Nx*image_size+start_depth)*1000, img_data)
axis('image')
set(gca,'FontSize',14)
colormap(gray(256))
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
axis([-50 50 20 105])

