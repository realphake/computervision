function noNoise = removeNoise(noise)
    noNoise = noise(noise(:,3) < 2,:);
end
