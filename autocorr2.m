% 2D autocorrelation calculations of the superconducting gap maps

function auto_corr_matrix = autocorr2(I)
% Using the Wiener-Khinchin theorem

I = double(I); % convert values to double
I = I - mean(I(:)); % substract the mean
I = I/sqrt(sum(I(:).^2)); % normalization
fft_I = fft2(I); % Fourier transform computation
auto_corr_matrix = real(fftshift(ifft2(fft_I.*conj(fft_I)))); % autocorrelation calculations