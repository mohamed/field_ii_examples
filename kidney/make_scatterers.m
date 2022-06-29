%  Make the scatterer file
%
%  Version 1.1, August 16, 2007, JAJ

N=1e6;
[phantom_positions, phantom_amplitudes] = human_kidney_phantom (N);

%  Save the data

save pht_data.mat phantom_positions phantom_amplitudes
