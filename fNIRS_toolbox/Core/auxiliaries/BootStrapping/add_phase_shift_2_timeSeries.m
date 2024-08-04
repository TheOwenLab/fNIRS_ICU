function[y_shifted] = add_phase_shift_2_timeSeries(y,limit)
% Code to generate Surrogate Data for Time Series 
% based on "phase-randomized" Fourier Transformation. 
% 
% Input: 
%   y - time-seris: nx1 vector in which n is the time points.
%   limit - limit = [a,b] vector with the shift interval
% Output: 
%   y_shifted: phase-randomized time series (Scrambled data).
%   

% Perform fft in the inputed data 
%(fast fourier transformation)
y_fft = fft(y);

% random phase between 0 and 2pi for each frequency
L = length(y);
shift = limit(1) + (limit(2)-limit(1)).*rand(L,1); 

% Option 1:
y_fft_shifted = y_fft.*exp(1i.*(shift));

y_shifted = real(ifft(y_fft_shifted));

end