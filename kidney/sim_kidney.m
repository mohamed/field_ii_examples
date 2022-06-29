%  Phased array B-mode scan of a human kidney
%
%  This script assumes that the field_init procedure has been called
%  Here the field simulation is performed and the data is stored
%  in rf-files; one for each rf-line done. The data must then
%  subsequently be processed to yield the image. The data for the
%  scatteres are read from the file pht_data.mat, so that the procedure
%  can be started again or run for a number of workstations.
%
%  Example by Joergen Arendt Jensen and Peter Munk, 
%  Version 1.1, April 1, 1998, JAJ.

%  Ver. 1.1: 1/4-98: Procedure xdc_focus_center inserted to use the new
%                    focusing scheme for the Field II program

%  Generate the transducer apertures for send and receive

f0=7e6;                  %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
lambda=c/f0;             %  Wavelength [m]
width=lambda/2;          %  Width of element
element_height=5/1000;   %  Height of element [m]
kerf=lambda/10;          %  Kerf [m]
focus=[0 0 90]/1000;     %  Fixed focal point [m]
N_elements=128;          %  Number of physical elements

%  Set the sampling frequency

set_sampling(fs);
set_field ('show_times', 5)

%  Generate aperture for emission

xmit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 5, focus);

%  Set the impulse response and excitation of the xmit aperture

impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (xmit_aperture, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (xmit_aperture, excitation);

%  Generate aperture for reception

receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 5, focus);

%  Set the impulse response for the receive aperture

xdc_impulse (receive_aperture, impulse_response);

%   Load the computer phantom

load pht_data

%  Set the different focal zones for reception

focal_zones=[5:1:150]'/1000;
Nf=max(size(focal_zones));
focus_times=(focal_zones-10/1000)/1540;
z_focus=60/1000;          %  Transmit focus

%  Set the apodization

apo=hanning(N_elements)';
xdc_apodization (xmit_aperture, 0, apo);
xdc_apodization (receive_aperture, 0, apo);

%   Do phased array imaging

no_lines=128;                  %  Number of lines in image
image_width=90/180*pi;         %  Size of image sector [rad]
dtheta=image_width/no_lines;   %  Increment for image

% Do imaging line by line

for i=1:128

  if ~exist(['rf_data/rf_ln',num2str(i),'.mat'])
    
    cmd=['save rf_data/rf_ln',num2str(i),'.mat i']
    eval(cmd)

  %   Set the focus for this direction

    theta= (i-1-no_lines/2)*dtheta;
    xdc_focus (xmit_aperture, 0, [z_focus*sin(theta) 0 z_focus*cos(theta)]);
    xdc_focus (receive_aperture, focus_times, [focal_zones*sin(theta) zeros(max(size(focal_zones)),1) focal_zones*cos(theta)]);
  
    %   Calculate the received response

    [rf_data, tstart]=calc_scat(xmit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);

    %  Store the result

    cmd=['save rf_data/rf_ln',num2str(i),'.mat rf_data tstart']
    eval(cmd)
    end

  end

%   Free space for apertures

xdc_free (xmit_aperture)
xdc_free (receive_aperture)

