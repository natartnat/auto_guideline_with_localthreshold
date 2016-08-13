%% turn to RGB
function [imt] = turnRGB(imt)
c = makecform('cmyk2srgb'); % Only need to be called once)
if size(imt,3) == 1 % Grayscale image
imt = cat(3, imt, imt, imt) ;
elseif size(imt,3) == 4 % CMYK image
imt(isnan(imt)) = 0;
imt = single(applycform(double(imt), c) .* 255); % Turn CMYK into RGB image
end