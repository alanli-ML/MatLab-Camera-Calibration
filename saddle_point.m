function [pt] = saddle_point(I)
% SADDLE_POINT Locate saddle point in an image patch.
%
%   [pt] = SADDLE_POINT(I) finds the subpixel center of a cross-junction in the 
%   image patch I, by blurring the patch, fitting a hyperbolic paraboloid to 
%   it, and then finding the critical point of that paraboloid.
%
%   Note that the location of 'p' is relative to (0.5, 0.5) at the upper left 
%   corner of the patch, i.e., the pixels are treated as covering an area of
%   one unit square.
%
%   Inputs:
%   -------
%    I  - mxn image patch (grayscale, double or integer class).
%
%   Outputs:
%   --------
%    pt  - 2x1 subpixel location of saddle point in I (x, y coords).
%
% References:
%
%   L. Lucchese and S. K. Mitra, "Using Saddle Points for Subpixel Feature
%   Detection in Camera Calibration Targets," in Proc. Asia-Pacific Conf.
%   Circuits and Systems (APCCAS'02), vol. 2, (Singapore), pp. 191-195,
%   Dec. 2002.


    %get size of input
    m = size(I,1);
    n = size(I,2);
   
    
    %apply gaussian blur
    kernel_size = 3;
    Ib = gaussian_blur(I,kernel_size,2);
    %crop padded areas from convolution; inaccurate data
    crop_pad = ceil(kernel_size/2);
    Ib = Ib(crop_pad+1:m-crop_pad,crop_pad+1:n-crop_pad);
    
    %get new size of input
    m = size(Ib,1);
    n = size(Ib,2);
    
    %additional filtering for non crossjunctions: a patch with a
    %crossjunction will have two white blobs and two black blobs. Search
    %along the edges of the image for these blobs.
%    imshow(Ib);
    %binarize image
    threshold = (max(Ib) + min(Ib))/2;
    I_bin = Ib;
    I_bin(I_bin<threshold) = 0;
    I_bin(I_bin>=threshold) = 1;
 %   imshow(I_bin*254);
    
    %concatenate sides of image into single continuous vector
    comb_sides = [(I_bin(1,:)),I_bin(2:end,n)',flip(I_bin(m,1:end-1)),flip(I_bin(2:end-1,1))'];
    
    %find a starting edge and rearrange the array; repeat last element to
    %check for change between first and last element
    for i=2:size(comb_sides,2)
         if comb_sides(i-1) ~= comb_sides(i)
             comb_sides = [comb_sides(i:end),comb_sides(1:i)];
         end
    end
    
    %moving along the sides of the image, it needs to change between black and white 4 times for
    %two black blobs to be present
    curr_pix = comb_sides(1);
    num_pix = 1;
    num_changes = 0;
    for i=2:size(comb_sides,2)
         if curr_pix ~= comb_sides(i)
             %if length of band is less 2, consider as noise
             if num_pix > 2
                 num_changes = num_changes + 1;
             end
             num_pix = 1;
         else
             num_pix = num_pix + 1;
         end
         curr_pix = comb_sides(i);
    end
    %no saddle point if less than four value changes
    if num_changes < 4
        pt = [-1;-1];
        return
    end
    %create matrix of x values
    X = zeros(m*n,6);
    %create vector of target y values
    Y = zeros(m*n,1);
    
    idx = 1;
    for x=1:n
        for y =1:m
            %assign each x vector for each point
            X(idx,:) = [x^2, x*y, y^2, x, y, 1];
            Y(idx) = Ib(y,x);
            idx = idx + 1;
        end
    end
    
    %calculate linear least squares directly using formula
    a = inv(X' * X)* X' * Y;
    
    %calulate saddle point using intersection of two lines
    A = [2*a(1) a(2); a(2) 2*a(3)];
    b = [ -a(4); -a(5)];
    pt = A\b;
    pt = pt + crop_pad;
  
%------------------
  
end