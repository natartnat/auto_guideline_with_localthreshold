function [localMaxImage,localMinImage] = bernsenTfilter(inputImage,winLength)
localMaxImage = colfilt(inputImage, [winLength winLength], 'sliding', @max);
localMinImage = colfilt(inputImage, [winLength winLength], 'sliding', @min);
end