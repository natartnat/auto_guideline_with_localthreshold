% testing rotate
function [] = CrossSectionConstruct()
clear all;close all;
% load datatube;
[filename,pathname,~] = uigetfile({'*.*';'*.png';'*.jpg';'*.tif'}, 'Get Image file');

tempname = filename;
% read path
fich1 = fullfile(pathname,tempname);[~,name,~] = fileparts(tempname);
im1 = imread(fich1);
imGray = rgb2gray(im1);
% save Load image
%    load A7Load.mat
ff = edge(imGray,'canny');
[y, x] = find(ff == 1);
edgesnap = [x(:),y(:)];
figure,imshow(im1);
edgesnap = sortedPoint(edgesnap);
% user mask
for i = 1:2
    k = 1;
    while k == 1
        [~, maskX,maskY] = roipoly;
        disp(['if proceed please click but if not please press a button']);
        k = waitforbuttonpress;
    end
    xy = [];
    for count = 1:size(maskX,1)-2
        %add data on field
        [x, y]=bresenham(maskX(count),maskY(count),maskX(count+1),maskY(count+1));
        xy = [xy;x,y];
    end
    % Snap line
    matching = knnclassify(xy,edgesnap,1:1:size(edgesnap,1));
    newSnap = edgesnap(matching,1:2);
    hold on,plot(edgesnap(matching,1),edgesnap(matching,2),'b*');
    distCross(i) = pdist([newSnap(1,:);newSnap(end,:)]);

    % Rotate line
    newSnap = rotateLine2(newSnap); % rotate for other half
    plot(abs(newSnap(:,1)),abs(newSnap(:,2)),'*g')

    Cross(i).Snapline = newSnap;

    if i == 1
        start = [Cross(1).Snapline(1,:);Cross(1).Snapline((end-1)/2,:)];
    elseif i==2
        stop  = [Cross(2).Snapline(1,:);Cross(2).Snapline((end-1)/2,:)];
    end
end
% profile on;
% load A3Cross;
%   load A6Cross
ws = [3 5 7 9 11 13 15 17 19 21 23 25 27 29 31];
%   load A8Cross
% fig = figure;
% load file
lowB = 0;
for j = 1:size(ws,2)
    A = calContrast(imGray,[ws(j) ws(j)],15);
    proposed = 'ourMethod'; % ourmethod
    % weak noise filter
    A=imfilter(A,fspecial('average',ws(j)+2),'replicate');
    % boundary and flood fill
    BW = imfill(A,'holes');
    [B,~,N,A] = bwboundaries(BW);
    
    % show result compare and error
    if size(B,1) == 1 ;
        % stop when has only 1 contour
        disp([name ' ' proposed ' W :' num2str(ws(j)) ...
            ' Nobj: ' num2str(length(B))  'Break'])
        break;
    elseif ws(j) >= 31
        disp([name ' ' proposed ' W :' num2str(ws(j)) ...
            ' Nobj: ' num2str(length(B)) ' nofound']);
        for g = 1:length(B)
            % assume more point is largest
            if length(B{g}) > lowB
                B{1} = B{g};
                lowB = length(B{g});
            end
        end
    end
end
% put B out of box
B = B{1};
B = [B(:,2),B(:,1)];
%% profile creation
% function XXXXX() Cross section replacement
%% start stop matching

for a = 1:2
    matching = knnclassify(start(a,:),B,1:1:size(B,1));
    start(a,:) = B(matching,:);
end
for a = 1:2
    matching = knnclassify(stop(a,:),B,1:1:size(B,1));
    stop(a,:) = B(matching,:);
end

%% ending touch
fig_out = figure;
[x_out,y_out,z_out] = rotationV3(B,BW,start,stop,Cross);
%% lighting & surface
% resort information
% for i = 1:size(x_out,1)
%     recon.x = [];
%     recon.y = [];
%     recon.z = [];
%     tempX = x_out{i};    tempY = y_out{i};    tempZ = z_out{i};
%     for c = 1:size(tempX,1)
%         recon.x = [recon.x ;tempX'];
%         recon.y = [recon.y ;tempY'];
%         recon.z = [recon.z ;tempZ'];
%     end
%     out(i).recon = recon;
% end
%create 3D 
% fig_out = figure;
% for i = 1:size(x_out,1)
%     temp = out(i).recon;
%     reconX = temp.x;reconY = temp.y;reconZ = temp.z;
%     surf(reconX,reconY,reconZ,'MeshStyle','column','FaceColor','b',...
%         'FaceLighting','gouraud','EdgeColor','b','EdgeAlpha',0.1);
%     hold on;
%     plot3((reconX(end,1:end))',(reconY(end,1:end))',(reconZ(end,1:end))','c');
%     fill3((reconX(end,1:end))',(reconY(end,1:end))',(reconZ(end,1:end))','b')
%     plot3((reconX(1,1:end))',(reconY(1,1:end))',(reconZ(1,1:end))','c');
%     fill3((reconX(1,1:end))',(reconY(1,1:end))',(reconZ(1,1:end))','b');
% end
savingPath = strcat(pathname,'resultCompare\model\');
Result_of = strcat('_result_','_ws',num2str(ws(j)));
saveas(fig_out,strcat(savingPath,name,Result_of,'.png'));
saveas(fig_out,strcat(savingPath,name,Result_of,'.fig'));
clf(fig_out);
profile viewer;
profile off;
end

function [localmask] = calContrast(image,window,contrast_threshold)

image = double(image);
localMax = maxfilt2(image, window);
localMin = minfilt2(image, window);
local_contrast = localMax - localMin;

localmask = local_contrast >= contrast_threshold;

end
