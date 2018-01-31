function [Ipts] = cross_junctions(I, boundPoly, Wpts)
% CROSS_JUNCTIONS Find cross-junctions in image with subpixel accuracy.
%
%   [Ipts] = CROSS_JUNCTION(I, boundPoly, Wpts) locates a series of cross- 
%   junction points on a planar calibration target, where the target is
%   bounded in the image by the specified 4-sided polygon. The number of
%   cross-junctions identified should be equal to the number of world points.
%
%   Note also that the world and image points must be in *correspondence*,
%   that is, the first world point should map to the first image point, etc.
%
%   Inputs:
%   -------
%    I          - Image (grayscale, double or integer class).
%    boundPoly  - 2x4 array defining bounding polygon for target (clockwise).
%    Wpts       - 3xn array of world points (in 3D, on calibration target).
%
%   Outputs:
%   --------
%    Ipts  - 2xn array of cross-junctions (x, y), relative to the upper left
%            corner of I.

    %Ib = gaussian_blur(I,5,3);
    Ib = I;
    Ipts = [];
    r = 7;
    [corners,strengths] = harris_corners(Ib);
    
    
    %apply non-maximum suppression
    idx = 1;
    dist_threshold = 15;
    while idx < size(corners,2)
       corner_max = idx;
       idx = idx + 1;
       if idx >= size(corners,2)
           break;
       end
       idx2 = idx;
       while idx2 <= size(corners,2)
           %remove corner if a stronger one is present within specified
           %threshold
           if(norm(corners(:,corner_max) - corners(:,idx2)) < dist_threshold)
               if strengths(corner_max) < strengths(idx2)
                   corners(:,corner_max) = corners(:,idx2);
                   strengths(corner_max) = strengths(idx2);
               end
                   corners(:,idx2) = [];
                   strengths(idx2) = [];
           else
               idx2 = idx2 + 1;
           end
           
       end
    end
    pts = corners;
    
    %get all corners within bounding polygon
    for p = 1:size(pts,2)
        if inpolygon(pts(1,p),pts(2,p),boundPoly(1,:),boundPoly(2,:)) 
            
            %create roi around corner
            corner = I(pts(2,p)-r:pts(2,p)+r,pts(1,p)-r:pts(1,p)+r);
            %find saddle point if present
            pt = saddle_point(corner);
            %remove saddle points too far from center of roi
            if pt(1) < r/3 || pt(1) > r + 2/3*r || pt(2) < r/3 || pt(2) > r + 2/3*r
                continue
            end
            
            %hold on
            %scatter(pt(1),pt(2))
            
            %get saddle point in image frame
            pt(1) = pt(1)+pts(1,p)-r-1;
            pt(2) = pt(2) + pts(2,p)-r-1;
            
            %imshow(uint8(Ib))
            %hold on
            %scatter(pt(1),pt(2))
            Ipts = cat(2,Ipts,pt);
        end
    end
    
    imshow(uint8(Ib))
    hold on 
    scatter(pts(1,:),pts(2,:))
    scatter(Ipts(1,:),Ipts(2,:))
    scatter([0,640],[0,480])
    
    
    %now to match points
    
    % For each corner of bounding polygon, find the closest saddle point
    % these are the corner points of the grid.
    corner_pts = zeros(2,4);
    for i=1:size(Ipts,2)
        for j = 1:4
            if norm(boundPoly(:,j) - corner_pts(:,j)) > norm(boundPoly(:,j) - Ipts(:,i))
                corner_pts(:,j) = Ipts(:,i); 
            end
        end
    end
        
    %check if first two corner points make short side or long side 
    side1 = corner_pts(:,1) - corner_pts(:,2);
    side2 = corner_pts(:,2) - corner_pts(:,3);
    side3 = corner_pts(:,3) - corner_pts(:,4);
    side4 = corner_pts(:,4) - corner_pts(:,1);
    
    %create a polgyon around each side to bound one row/column of points
    side_bound1 = [corner_pts(:,1) + 0.05*side4, corner_pts(:,1) - 0.05*side4, ...
        corner_pts(:,2) - 0.05*side2, corner_pts(:,2) + 0.05*side2];
    
    side_bound2 = [corner_pts(:,2) + 0.05*side1, corner_pts(:,2) - 0.05*side1, ...
        corner_pts(:,3) - 0.05*side3, corner_pts(:,3) + 0.05*side3];
    
    %get number of points contained in each polygon
    side1_numPoints = sum(inpolygon(Ipts(1,:),Ipts(2,:),side_bound1(1,:),side_bound1(2,:)));
    side2_numPoints = sum(inpolygon(Ipts(1,:),Ipts(2,:),side_bound2(1,:),side_bound2(2,:)));
    
    if side1_numPoints < side2_numPoints
        %swap to make sure long side comes first
        corner_pts = [corner_pts,corner_pts(:,1)];
        corner_pts(:,1) = [];
    end
   
    
    
    %assuming world points are in 6x8 grid with identical Z values:
    corner_pts_world = [Wpts(1:2,1),Wpts(1:2,8),Wpts(1:2,48),Wpts(1:2,41)];
    
    %find homography between corner saddle points and world points    
    H = dlt_homography(corner_pts_world,corner_pts);
    Wpts_trans = H * [Wpts(1:2,:);ones(1,size(Wpts,2))];
    Wpts_trans = Wpts_trans(1:2,:) ./ Wpts_trans(3,:);
    %find nearest neighbor for each world point
    for i=1:size(Wpts_trans,2)
        for j=i:size(Ipts,2)
            %if current point is closer to world, swap
            if norm(Wpts_trans(:,i)-Ipts(:,i)) > norm(Wpts_trans(:,i)-Ipts(:,j))
                 temp = Ipts(:,i);
                 Ipts(:,i) = Ipts(:,j);
                 Ipts(:,j) = temp;
            end
        end
    end
    
    %remove any extra unmatched points
    Ipts = Ipts(:,1:size(Wpts,2));
%------------------
  
end


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

%--- FILL ME IN ---

% Code goes here...
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

