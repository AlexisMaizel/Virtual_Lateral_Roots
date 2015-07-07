function real = convertToReal( value )

% invert the log
real = exp( value );
% multiply by length of pixel in µm
real = real * 0.3225;
% divide by temporal difference in minutes
real = real / 63.5;
% and take again the log
real = log10( real );
