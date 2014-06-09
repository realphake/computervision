function [F,matches, f1,f2] = FME( impath1, impath2, method )
    run('vlfeat-0.9.18/toolbox/vl_setup');
    im1 = imresize(imread(impath1), [600 NaN]);
    im2 = imresize(imread(impath2), [600 NaN]);
    im1 = imcrop(im1,[250,120,330,330]);
    im2 = imcrop(im2,[250,120,330,330]);
    [f1, d1] = vl_sift(single(rgb2gray(im1)));
    [f2, d2] = vl_sift(single(rgb2gray(im2)));
    [matches, ~] = vl_ubcmatch(d1, d2);
    if ( strcmp(method, 'ransac') )
        [F,bestSetOfInliers] = FMEransac( f1, f2, matches, 0.1 );
    elseif ( strcmp(method, 'normalized') )
        F = FMEnorm( f1, f2, matches );
    elseif ( strcmp(method, 'regular') )
        F = FMEregular( f1, f2, matches );
    else
        F = eye(3);
    end
    matches = bestSetOfInliers;
    
end

function [F,bestSetOfInliers] = FMEransac( f1, f2, matches, threshold )
    bestSetOfInliers = [];
    for iterations = 1:1000
        subsample = randomcolumns(matches, 8);
        Ftemp = FMEnorm( f1, f2, subsample );
        theseInliers = [];
        for match = matches
            p = f1(:,match(1));
            pn = f2(:,match(2));
            p = [p(1:2);1];
            pn = [pn(1:2);1];
            numerator = (pn'*Ftemp*p)^2;
            Fp = (Ftemp*p);
            Fpn = (Ftemp'*pn);
            denominator = Fp(1)^2 + Fp(2)^2 + Fpn(1)^2 + Fpn(2)^2;
            distance = numerator / denominator;
            if (distance < threshold)
                theseInliers = [theseInliers,match];
            end
        end
        if (length( theseInliers ) > length( bestSetOfInliers ))
            bestSetOfInliers = theseInliers;
        end
    end
    F = FMEnorm( f1, f2, bestSetOfInliers );
end

function picks = randomcolumns(columns, number)
    picks = zeros(2,number);
    for i = 1:number
        picks(:,i) = columns( :, random('unid', length(columns)) );
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
    Fhat = FMEregular( f1hat, f2hat, matches );
    F = T2' * Fhat * T1;
    
end
