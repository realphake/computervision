function [ output_args ] = FME( impath1, impath2 )
    im1 = imread(impath1);
    im2 = imread(impath2);
    [f1, d1] = vl_sift(single(im1));
    [f2, d2] = vl_sift(single(im2));
    [matches, scores] = vl_ubcmatch(d1, d2);

end
