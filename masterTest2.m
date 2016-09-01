function [errorA, errorOur] = masterTest2()
clear all;close all;clc;
profile on;
%      [filename,pathname,~] = uigetfile({'*.*';'*.png';'*.jpg';'*.tif'}, 'Get Image file','MultiSelect','on');
%      [filenameAns,pathnameAns,~] = uigetfile({'*.*';'*.png';'*.jpg';'*.tif'}, 'Get Answer file','MultiSelect','on');
%      load loadimage.mat;
load newImageNo2.mat
T = 15; k = 0.06 ;kSau = 0.06;
ws = [3 5 7 9 11 13 15 17 19 21 23 25 27 29 31];
if iscell(filename)
    loop = size(filename,2);
else
    loop = 1;
end
%     load toothpasteCompare
fig = figure;
for type=1:6
    for i = 1:1:loop
        % load file
        try
            tempname = filename{i};
            tempAns = filenameAns{i};
        catch
            %if only 1 files
            tempname = filename;
            tempAns = filenameAns;
        end
        % read path
        fich1 = fullfile(pathname,tempname);[~,name,~] = fileparts(tempname);
        fichAns = fullfile(pathnameAns,tempAns);
        % save path
        savingPath = strcat(pathname,'resultCompare\ourMethod_newImage_01092016_noisereduct\');
        % convert to gray
        try
            im1 = imread(fich1);
            imGray = rgb2gray(im1);
            imAns = imread(fichAns);
            imGrayAns = rgb2gray(imAns);
        catch
            im1 = imread(fich1);
            imGray = mat2gray(im1);
            % if object is not white
            if sum(sum(imGray))/(size(im1,1)*size(im1,2)) < 1
                imGray = imcomplement(imGray);
            end
            imAns = imread(fichAns);
            imGrayAns = mat2gray(imAns);
        end
        % Loop for size
        for j = 1:size(ws,2)
            switch type
                case 1 ,A = calContrast(imGray,[ws(j) ws(j)],T);
                    proposed = 'ourMethod'; % ourmethod
                case 2 ,A = imcomplement(im2bw(adaptivethreshold(imGray,ws(j),k,0)));
                    proposed = 'mean';  % mean
                case 3 ,A = imcomplement(im2bw(adaptivethreshold(imGray,ws(j),k,1)));
                    proposed = 'median';% median
                case 4 ,A = imcomplement((sauvola(imGray,[ws(j) ws(j)],kSau)));
                    proposed = 'sauvola';% sauvola
                case 5 ,A = imcomplement(im2bw(NEWLOCAL(imGray,[ws(j) ws(j)],k)));
                    proposed = 'Singh';% new local T.R Singh’s
                case 6 ,A = imcomplement(im2bw(bernsenE2(imGray,[ws(j) ws(j)],k)));
                    proposed = 'bernsen';% bernsen
            end
            % naming path
            Result_of = strcat('_result_',proposed,'_localContrast_T',num2str(T),'_ws',num2str(ws(j)));
            % size of image
            [M ,N] = size(A);
            BeReduce = [];
            % boundary and flood fill
            [BW,B,~,~] = fillMask(A);
            [BWAns,~,~,~] = fillMask(imGrayAns);
            % error compare
            error = (sum(sum(abs( im2bw(BW)-im2bw(BWAns) ))) / (M*N)) * 100;
            % show comparison
            subplot(1,2,1),
            imshow(BW);title(['windowsize : ' num2str(ws(j)) 'x' num2str(ws(j)) ' error : ' num2str(error)]);
            % shutdown warning
            try
                w = warning('query','last');
                id = w.identifier;
                warning('off',id);
            catch
            end
            % [underConstruction]reduce noise
            % BeReduce = length(B);
            % [B] = areaFilter(B,M*N,0.10);
            
            % show line for our method
            if length(B) == 1
                boundary = B{1};
                hold on,plot(boundary(:,2), boundary(:,1),'LineWidth',2);
            end
            
            % show in image
            subplot(1,2,2),imshow(abs(BW-im2bw(BWAns)));title(['error : ' num2str(error)]);
            
            % save image
            saveas(fig,strcat(savingPath,name,Result_of,'.png'));
            % clear fig
            clf(fig);
            
            % show result compare and error
            if size(B,1) == 1 ;
                % stop when has only 1 contour
                disp([name ' ' proposed ' err: ' num2str(error) ' W :' num2str(ws(j)) ...
                    ' Nobj: ' num2str(length(B)) ' Bobj: ' num2str(BeReduce) ' Break'])
                break;
            elseif ws(j) >= 31
                disp([name ' ' proposed ' err: ' num2str(error) ' W :' num2str(ws(j)) ...
                    ' Nobj: ' num2str(length(B)) ' Bobj: ' num2str(BeReduce) ' nofound'])
            else
%                 disp([name ' ' proposed ' err: ' num2str(error) ' W :' num2str(ws(j)) ...
%                     ' Nobj: ' num2str(length(B)) ' Bobj: ' num2str(BeReduce) ])
            end
        end
    end
end
end
function [BW,B,L,N] = fillMask(image)
BW = imfill(image,'holes');
[B,L,N] = bwboundaries(BW);
end
function [localmask] = calContrast(image,window,contrast_threshold)

image = double(image);
localMax = maxfilt2(image, window);
localMin = minfilt2(image, window);
local_contrast = localMax - localMin;

localmask = local_contrast >= contrast_threshold;

end

function [Breturn] = areaFilter(Boundary,imageSize,condition)
% return reduction of noise area
Breturn = {};num = 1;
for i = 1:length(Boundary)
    matrix = cell2mat(Boundary(i));
    area = polyarea(matrix(:,1),matrix(:,2));
    if area < condition * imageSize
        continue;
    else
        Breturn(num) = Boundary(i);
        num = num + 1;
    end
end

end
%% dump code zone
function []  = plotCompare(im1,A,imGrayAns,j)
% unused
subplot(2,3,1),imshow(im1);
title(['original image']);
[M ,N] = size(A);
[BW,B,~,~] = fillMask(A);
[BWAns,~,~,~] = fillMask(imGrayAns);
% fill mask
subplot(2,3,2),imshow(BW);
title(['mask : ' num2str(ws(j))]);
% subtract mask
subplot(2,3,3),imshow(im2bw(BWAns));
title(['manual mask']);
subplot(2,3,4),imshow(BW-im2bw(BWAns));
error = abs(sum(sum(BW-im2bw(BWAns)))/(M*N))/100;
title(['error' num2str(error)]);
subplot(2,3,5),imshow(im1);
title(['window : ' num2str(ws(j)) 'number contour']);
end
function [Breturn] = noiseReduction(M,N,boundary,B)
%objArea = (max(boundary(:,1))-min(boundary(:,1)))*(max(boundary(:,2))-min(boundary(:,2)));
%['objArea : ' num2str(objArea/(M*N))]
%if objArea > 0.01*M*N %|| length(B) == 1;
%ans = B(i);
%break;
%end
end
