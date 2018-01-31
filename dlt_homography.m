function [H, A] = dlt_homography(I1pts, I2pts)
% dlt_homography Perspective Homography between two images.
%
%   Given 4 points from 2 separate images, compute the perspective homograpy
%   (warp) between these points using the DLT algorithm.
%
%   Inputs:
%   -------
%    I1pts  - 2x4 array of points from Image 1 (each column is x, y).
%    I2pts  - 2x4 array of points from Image 2 (1-to-1 correspondence).
%
%   Outputs:
%   --------
%    H  - 3x3 perspective homography (matrix map) between image coordinates.
%    A  - 8x9 DLT matrix used to determine homography.

  A = zeros(8,9);
  %fill in matrix
  for i = 1:4
     A(2*i-1,:) =  [-I1pts(1,i),-I1pts(2,i),-1,0,0,0,I2pts(1,i)*I1pts(1,i),I2pts(1,i)*I1pts(2,i),I2pts(1,i)];
     A(2*i,:) = [0,0,0,-I1pts(1,i),-I1pts(2,i),-1,I2pts(2,i)*I1pts(1,i),I2pts(2,i)*I1pts(2,i),I2pts(2,i)];
  end
  %find null space
  n = null(A);
  H = reshape(n,[3,3]);
  H = H';
  %disp(H)
%------------------
  
end
