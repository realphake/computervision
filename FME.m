% This is the wrapper function. It calls one of the three functions below,
% which actually call each other hierarchically. Robert Martin might have
% some words to say about that, but it's simple enough.
function [F,matches, f1,f2] = FME( impath1, impath2, method )
	% Installing this in the code for convenience.
    run('vlfeat-0.9.18/toolbox/vl_setup'); 
    
	tic
	% This is kind of specific to the image sizes that you're using. These
	% numbers are meant for the house dataset. 
    % The bear has to be cropped, too.
    im1 = imresize(imread(impath1), [300 NaN]); 
    im2 = imresize(imread(impath2), [300 NaN]); % Loaded and resized.
    [f1, d1] = vl_sift(single(im1));
    [f2, d2] = vl_sift(single(im2)); % SIFT features are extracted.
	[matches, ~] = vl_ubcmatch(d1, d2); % Matches have been found.
	
	% Here the method variable is used to determine what function will be 
    % called. Each of these functions calls the one below it, too. 
    if ( strcmp(method, 'ransac') )
		% Note this method changes which matched SIFT features will be 
        % returned. Specifically, it will only return the inliers.
        [F,matches] = FMEransac( f1, f2, matches, 0.1 );
    elseif ( strcmp(method, 'normalized') )
        F = FMEnorm( f1, f2, matches );
    elseif ( strcmp(method, 'regular') )
        F = FMEregular( f1, f2, matches );
    else
        F = eye(3); % We should not get here, of course.
    end
    toc
end

% This subfunction determines the set of inliers and then calls the
% normalized variant of the function on that set. Also returns the
% set of inliers.
function [F,bestSetOfInliers] = FMEransac( f1, f2, matches, threshold )
    bestSetOfInliers = [];
    for iterations = 1:100 % Takes about 6-9 seconds.
        % Take 8 random matches, calculate a fundamental matrix from them
        % using the normalized method.
        subsample = randomcolumns(matches, 8);
        Ftemp = FMEnorm( f1, f2, subsample );
        theseInliers = [];
        % Find how many total points behave according to that matrix. 
        % Save the largest number of inliers found.
        for match = matches
            
            % p and pn are the coordinates of this match in the first and
            % second frames, respectively.
            p = f1(:,match(1));
            pn = f2(:,match(2));
            p = [p(1:2);1];
            pn = [pn(1:2);1];
            
            % This is the calculation for the distance between a point and
            % where it should be according to the fundamental matrix.
            numerator = (pn'*Ftemp*p)^2;
            Fp = (Ftemp*p);
            Fpn = (Ftemp'*pn);
            denominator = Fp(1)^2 + Fp(2)^2 + Fpn(1)^2 + Fpn(2)^2;
            distance = numerator / denominator;
            
            % If that distance is small enough, it's an inlier.
            if (distance < threshold)
                theseInliers = [theseInliers,match];
            end % end if (distance < threshold)
            
        end % end for match = matches
        
        if (length( theseInliers ) > length( bestSetOfInliers ))
            % If this set of inliers is larger than what we've found so
            % far, we're going to replace our best candidate.
            bestSetOfInliers = theseInliers;
        end
    end
    % Finally calculate the fundamental matrix one last time.
    F = FMEnorm( f1, f2, bestSetOfInliers );
end

% Pay no mind. This is just to find 8 random picks from the matches.
function picks = randomcolumns(columns, number)
    picks = zeros(2,number);
    for i = 1:number
        picks(:,i) = columns( :, random('unid', length(columns)) );
    end
end

% And here's the variant of the algorithm that only normalizes the
% features it uses. It's a buncha math.
function F = FMEnorm( f1, f2, matches )
    
    f1 = [f1(1:2,:); ones(1,length(f1))];
    mx = sum( f1(1,:) ) / length(f1); % Coulda just used mean().
    my = sum( f1(2,:) ) / length(f1);
    % d is the average distance of features to their average. I think?
    d = sum( sqrt( (f1(1,:)-mx).^2 + (f1(2,:)-my).^2 ) ) / length(f1);
    % And here I'm really forgetting things. I did this a couple weeks ago.
    T1 = [sqrt(2)/d, 0,         -mx*sqrt(2)/d;
          0,         sqrt(2)/d, -mx*sqrt(2)/d;
          0,         0,         1             ];
    
    % Anyway, do the same thing to the features of the second image.
    f2 = [f2(1:2,:); ones(1,length(f2))];
    mx = sum( f2(1,:) ) / length(f2);
    my = sum( f2(2,:) ) / length(f2);
    d = sum( sqrt( (f2(1,:)-mx).^2 + (f2(2,:)-my).^2 ) ) / length(f2);
    T2 = [sqrt(2)/d, 0,         -mx*sqrt(2)/d;
          0,         sqrt(2)/d, -mx*sqrt(2)/d;
          0,         0,         1             ];
    
    % So f-hat are the normalized features. They are found by
    % multiplication with the T-matrices.
    f1hat = T1*f1;
    f2hat = T2*f2;
    % And then we use the normal FME calculation from below, but on the
    % normalized features.
    Fhat = FMEregular( f1hat, f2hat, matches );
    % Un-normalize the result and done.
    F = T2' * Fhat * T1;
    
end

% And this, finally, is the actual calculation you need to do to get a
% fundamental matrix. Of course, the featuresets and matches are probably
% heavily edited by RANSAC and normalization, but it can be used without
% those things too.
function F = FMEregular( f1, f2, matches )

    % A is going to be a matrix with 9 columns and as many rows as there
    % are matches. The content of this matrix is not going to be any
    % clearer from a comment than if you just read the code.
    A = zeros(length(matches),9);
    for m = 1:length(matches)
        x1 = f1(1,matches(1,m)); 
        y1 = f1(2,matches(1,m)); % The x and y coordinates in image 1.
        x2 = f2(1,matches(2,m));
        y2 = f2(2,matches(2,m)); % The x and y coordinates in image 2.
        A(m,:) = [x1*x2, x1*y2, x1, y1*x2, y1*y2, y1, x2, y2, 1];
    end
    
    % Now A is used to calculate F, the fundamental matrix. From the third
    % component of a single value decomposition of A, we take the column
    % that has the smallest single element in it. That column,
    % restructured, is F.
    [~, ~, V] = svd(A);
    [~, column] = min(min(abs(V)));
    F = zeros(3,3);
    F(:,1) = V(1:3,column);
    F(:,2) = V(4:6,column);
    F(:,3) = V(7:9,column);
    
    % Here is a little refinement. The second element of another single
    % value decomposition is changed by setting the smallest element in its
    % diagonal to zero. Then F is reconstructed.
    [Uf, Df, Vf] = svd(F);
    diagonal = diag(Df);
    [~, columnf] = min(abs(diagonal));
    diagonal(columnf) = 0;
    newDf = diag(diagonal);
    F = Uf * newDf * Vf';
    
end
