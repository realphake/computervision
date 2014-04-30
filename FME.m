function F = FME( impath1, impath2, norm )
    im1 = imread(impath1);
    im2 = imread(impath2);
    [f1, d1] = vl_sift(single(im1));
    [f2, d2] = vl_sift(single(im2));
    [matches, scores] = vl_ubcmatch(d1, d2);
    
    if ( norm )
        F = FMEnorm( f1, f2, matches );
    else
        F = FMEregular( f1, f2, matches );
    end
end

function F = FMEregular( f1, f2, matches )
    A = zeros(length(matches),9);
    for m = 1:length(matches)
        x1 = f1(1,matches(1,m));
        y1 = f1(2,matches(1,m));
        x1p = f2(1,matches(2,m));
        y1p = f2(2,matches(2,m));
        A(m,:) = [x1*x1p, x1*y1p, x1, y1*x1p, y1*y1p, y1, x1p, y1p, 1];
    end
    [~, ~, V] = svd(A);
    [~, column] = min(min(abs(V)));
    F = zeros(3,3);
    F(:,1) = V(1:3,column);
    F(:,2) = V(4:6,column);
    F(:,3) = V(7:9,column);
    
    [Uf, Df, Vf] = svd(F);
    diagonal = diag(Df);
    [~, columnf] = min(abs(diagonal));
    diagonal(columnf) = 0;
    newDf = diag(diagonal);
    
    F = Uf * newDf * Vf';
    
end

function F = FMEnorm( f1, f2, matches )
    f1 = [f1(1:2,:); ones(1,length(f1))];
    mx = sum( f1(1,:) ) / length(f1);
    my = sum( f1(2,:) ) / length(f1);
    d = sum( sqrt( (f1(1,:)-mx).^2 + (f1(2,:)-my).^2 ) ) / length(f1);
    T1 = [sqrt(2)/d, 0,         -mx*sqrt(2)/d;
          0,         sqrt(2)/d, -mx*sqrt(2)/d;
          0,         0,         1             ];
    f2 = [f2(1:2,:); ones(1,length(f2))];
    mx = sum( f2(1,:) ) / length(f2);
    my = sum( f2(2,:) ) / length(f2);
    d = sum( sqrt( (f2(1,:)-mx).^2 + (f2(2,:)-my).^2 ) ) / length(f2);
    T2 = [sqrt(2)/d, 0,         -mx*sqrt(2)/d;
          0,         sqrt(2)/d, -mx*sqrt(2)/d;
          0,         0,         1             ];
    
    f1hat = T1*f1;
    f2hat = T2*f2;
    A = zeros(length(matches),9);
    for m = 1:length(matches)
        x1 = f1hat(1,matches(1,m));
        y1 = f1hat(2,matches(1,m));
        x1p = f2hat(1,matches(2,m));
        y1p = f2hat(2,matches(2,m));
        A(m,:) = [x1*x1p, x1*y1p, x1, y1*x1p, y1*y1p, y1, x1p, y1p, 1];
    end
    [~, ~, V] = svd(A);
    [~, column] = min(min(abs(V)));
    Fhat = zeros(3,3);
    Fhat(:,1) = V(1:3,column);
    Fhat(:,2) = V(4:6,column);
    Fhat(:,3) = V(7:9,column);
    
    [Uf, Df, Vf] = svd(Fhat);
    diagonal = diag(Df);
    [~, columnf] = min(abs(diagonal));
    diagonal(columnf) = 0;
    newDf = diag(diagonal);
    
    Fhat = Uf * newDf * Vf';
    F = T2' * Fhat * T1;
    
end


