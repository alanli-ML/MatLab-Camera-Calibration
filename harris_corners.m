function [corners,strength] = harris_corners(I)
    k = 0.04;
    Threshold = 350000000;
    sigma = 1;
    halfwid = sigma * 3;
    corners = [];
    strength = [];
    [xx, yy] = meshgrid(-halfwid:halfwid, -halfwid:halfwid);

    Gxy = exp(-(xx .^ 2 + yy .^ 2) / (2 * sigma ^ 2));

    Gx = xx .* exp(-(xx .^ 2 + yy .^ 2) / (2 * sigma ^ 2));
    Gy = yy .* exp(-(xx .^ 2 + yy .^ 2) / (2 * sigma ^ 2));


    numOfRows = size(I, 1);
    numOfColumns = size(I, 2);

    % 1) Compute x and y derivatives of image
    Ix = conv2(I,Gx,'same');
    Iy = conv2(I,Gy,'same');

    size(Ix)

    % 2) Compute products of derivatives at every pixel
    Ix2 = Ix .^ 2;
    Iy2 = Iy .^ 2;
    Ixy = Ix .* Iy;

    % 3)Compute the sums of the products of derivatives at each pixel
    Sx2 = conv2(Ix2,Gxy,'same');
    Sy2 = conv2(Iy2,Gxy,'same');
    Sxy = conv2(Ixy,Gxy,'same');

    im = zeros(numOfRows, numOfColumns);
    for x=1:numOfRows
       for y=1:numOfColumns

           % 4) Define at each pixel(x, y) the matrix H
           H = [Sx2(x, y) Sxy(x, y); Sxy(x, y) Sy2(x, y)];

           % 5) Compute the response of the detector at each pixel
           R = det(H) - k * (trace(H) ^ 2);

           % 6) Threshold on value of R
           if (R > Threshold)
               corners = cat(2,corners,[y;x]);
               strength = cat(2,strength,R);
              im(x, y) = R; 
           end
       end
    end

end