function real = convertToReal( value )

% invert the log
real = exp( value );
% multiply by length of pixel in µm
real = real * 0.3225;
% divide by temporal difference in minutes
real = real / 60; % corresponds to average time difference of 12.09 time steps
% and take again the log
real = log10( real );
