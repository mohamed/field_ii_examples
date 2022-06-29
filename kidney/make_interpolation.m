%  Function for making the interpolation of an ultrasound image.
%  The routine make_tables must have been called previously.
%
%  Input parameters: 
%
%         envelope_data - The envelope detected and log-compressed data as
%                         an integer array as 8 bits values
%
%  Output:  img_data    - The final image as 8 bits values
%
%  Calling: img_data = img_datamake_interpolation (envelope_data); 
%
%  Version 1.0, 14/2-1999, JAJ

function img_data = make_interpolation (envelope_data)

%  Call the appropriate function


img_data = fast_int (2, envelope_data);
	    
