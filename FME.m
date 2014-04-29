function output = FME( impath1, impath2 )
    im1 = imread(impath1);
    im2 = imread(impath2);
    [f1, d1] = vl_sift(single(im1));
    [f2, d2] = vl_sift(single(im2));
    [matches, scores] = vl_ubcmatch(d1, d2);
    
    A = zeros(length(matches),9);
    for m = 1:length(matches)
        x1 = f1(1,matches(1,m));
        y1 = f1(2,matches(1,m));
        x1p = f2(1,matches(2,m));
        y1p = f2(2,matches(2,m));
        A(m,:) = [x1*x1p, x1*y1p, x1, y1*x1p, y1*y1p, y1, x1p, y1p, 1];
    end
    [U, D, V] = svd(A);
    [smallest, column] = min(min(abs(V)));
    F = zeros(3,3);
    F(:,1) = V(1:3,column);
    F(:,2) = V(4:6,column);
    F(:,3) = V(7:9,column);
    
end
