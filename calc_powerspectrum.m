function [freq, x_ps, y_ps, x_0, y_0] = GetPowerSpectrum(t,x,y)

% Returns the power spectrum of the times series x and y and the equivalent
% frequency vector. t,x,y should be column vectors
%
% Note: GetPowerSpectrum was updated to rescale 2pi*old value on 090902, to
% give correct normalisation of single sided signals. Note that if only one
% side of a single sided spectrum is kept, the powers need to be doubled subsequent to this function.

% If less than 3 inputs, fix with zeros
if nargin < 3
    y = zeros( size( x ) ) ;
end

% If x, y and t are not all the same size return error
if( size( t , 1 ) ~= size( x , 1 ) | size( t , 1 ) ~= size( y , 1 ) )
    error( 'Vectors t,x,y are not of equal size' )
end

n     = length( x );
fNyq  = n / max( t ) / 2 ;                              % Nyquist frequency
x_fft = fft( x ); y_fft = fft( y );                     % Fourier transform data
x_0   = fftshift( x_fft ); y_0 = fftshift( y_fft );     % Shift values to centre on zero freq
x_ps  = x_0 .* conj( x_0 ); y_ps = y_0 .* conj( y_0 );  % Power spectra
freq  = ( [ 1 : n ] / max( t ) );                       % 0-centered frequency range

ind   = find( freq <= fNyq );                           % only to the Nyquist f
freq  = rot90( freq( ind ) );                           % rotate to keap everthing n x 1
x_ps  = x_ps( ind ); y_ps = y_ps( ind );


end

