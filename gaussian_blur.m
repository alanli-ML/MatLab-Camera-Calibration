function [Ib] = gaussian_blur(I, wndSize, sigma)
% GAUSSIAN_BLUR Smooth image with symmetric Gaussian filter.
%
%   [Ib] = GAUSSIAN_BLUR(I, wndSize, sigma) produces a filtered image Ib 
%   from I using a square Gaussian kernel with window size wndSize.
%
%   Inputs:
%   -------
%    I        - mxn intensity image.
%    wndSize  - Kernel window size (square, odd number of pixels).
%    sigma    - Standard deviation of Gaussian (pixels, symmetric).
%
%   Outputs:
%   --------
%    Ib  - mxn filtered output image, of same size and class as I.



    %create kernel
    kernel = zeros(wndSize);
    r = (wndSize - 1)/2;
    for i=-r:r
       for j=-r:r
           %get value of gaussian at x,y
           kernel(i+r+1,j+r+1) = gaussian2d(i,j,sigma);
       end
    end
    %normalize kernel to 1
    kernel = kernel./sum(reshape(kernel,wndSize^2,1));

    %------------------
    
    %apply convolution
    Ib = conv2(I,kernel, 'same');
    
    %convert data back to uint8 if input was uint8
    if isa(I,'uint8')
        Ib = uint8(Ib);
    end
end

%returns value of 2d gaussian with std sigma at x,y 
function g = gaussian2d(x,y,sigma)
    g = 1/(sqrt(2*pi)*sigma) * exp(-1/(2*sigma^2) * (x^2 + y^2));
end

