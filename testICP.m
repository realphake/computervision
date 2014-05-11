function [ ] = testICP( )
%TESTICP Summary of this function goes here
%   Detailed explanation goes here
    foo0 = readPcd('data/0000000000.pcd');
    foo1 = readPcd('data/0000000000.pcd');
    denoised0 = removeNoise(foo0);
    denoised1 = removeNoise(foo1);
    in0 = denoised0(:,1:3);
    in1 = denoised1(:,1:3);
    [output, R, T] = ICP(in0,in1,10,'random');
    R
    T
    
end

