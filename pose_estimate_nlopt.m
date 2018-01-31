function [E] = pose_estimate_nlopt(Eg, Ipts, Wpts)
%  POSE_ESTIMATE_NLOPT Estimate camera pose from 2D-3D correspondences via NLS.
%
%   [E] = POSE_ESTIMATE_NLOPT(Eg, Ipts, Wpts) performs a nonlinear least squares 
%   optimization procedure to determine the best estimate of the camera pose in 
%   the calibration target frame, given 2D-3D point correspondences.
%
%   Inputs:
%   -------
%    Eg    - 4x4 homogenous pose matrix, initial guess for camera pose.
%    Ipts  - 2xn array of cross-junction points (with subpixel accuracy).
%    Wpts  - 3xn array of world points (one-to-one correspondence with image).
%
%   Outputs:
%   --------
%    E  - 4x4 homogenous pose matrix, estimate of camera pose in target frame.

    N = size(Wpts,2);
    E = Eg;
    num_it = 100;
    lambda = 1;
    %load camera intrinsics
    intrinsic = [564.9 0 337.3; 0 564.3 226.5; 0 0 1];
    for i=1:num_it
        %create Hessian
        A = zeros(6);
        b = zeros(6,1);
        
        %convert camera pose to 3 rotation and 3 translation params
        S = [E(1:3,4);rpy_from_dcm(E(1:3,1:3))];
        
        %transform world points to image plane
        C = inv(E);
        y = C*[Wpts;ones(1,N)];
        y = [y(1,:)./y(3,:); y(2,:)./y(3,:); ones(1,N)];
        y = intrinsic* y;
        t = [Ipts(1,:); Ipts(2,:)];
        
        %get error between transformed points and detected cross junctions
        err = t-y(1:2,:);
        %sum(norm(err,1),2);
        %disp(norm(err));
        
        %calculate hessian over all points
        for j = 1:N
            cpt = Wpts(:,j);
            J = find_jacobian(intrinsic,E,cpt);
            A = A + J'*J;
            b = b + J'*err(:,j);
        end
        
        %solve for update p
        p = inv(A+lambda*diag(A))*b;
        %update previous estimation S
        S = S + p;
        %convert back to pose matrix
        E(1:3,1:3) = dcm_from_rpy(S(4:6));
        E(1:3,4) = S(1:3);
        %decrease lambda
        lambda = lambda/1.1;
    end
    
%------------------

end
