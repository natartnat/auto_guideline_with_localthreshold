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
    xy = []; clear s_in;
    while ~exist('s_in')|| s_in == 'e'
        [~, maskX,maskY] = roipoly;
        disp([' C:Circle , R:Rectangle, L:Line, E:Cancel']);
        s_in = input('request answer:','s' );
    end
    hold on;
    switch s_in
        case 'c',
            r = norm([maskX(1),maskY(1)] - [maskX(end-1),maskY(end-1)])/2 ;
            x = floor((maskX(1)+maskX(end-1))/2);
            y = floor((maskY(1)+maskY(end-1))/2);
            th = 0:pi/50:2*pi;
            xunit = r * cos(th) + x;
            yunit = r * sin(th) + y;
            xy = [xunit(:),yunit(:)];
        case 'r',
            for count = 1:size(maskX,1)-2
                %add data on field
                [x, y]=bresenham(maskX(count),maskY(count),maskX(count+1),maskY(count+1));
                xy = [xy;x,y];
            end
            xy = rotateLine2(xy);
            hold on;plot(xy(:,1),xy(:,2));
        case 'l',
            for count = 1:size(maskX,1)-2
                %add data on field
                [x, y]=bresenham(maskX(count),maskY(count),maskX(count+1),maskY(count+1));
                xy = [xy;x,y];
            end
        otherwise
            % user's defined
            for count = 1:size(maskX,1)-2
                %add data on field
                [x, y]=bresenham(maskX(count),maskY(count),maskX(count+1),maskY(count+1));
                xy = [xy;x,y];
                
            end
    end
       
    % Snap line
    Cross(i).Snapline = xy;
    
    if i == 1
        start = [maskX(1),maskY(1);maskX(end-1),maskY(end-1)];
    elseif i==2
        stop  = [maskX(1),maskY(1);maskX(end-1),maskY(end-1)];
    end
end
hold on; plot(xy(:,1),xy(:,2));
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
% Cross section replacement
xyz = [Cross(1).Snapline(:,1),Cross(1).Snapline(:,2),zeros(size(Cross(1).Snapline,1),1)];
xmean = mean(xyz(:,1)); ymean = mean(xyz(:,2));
xyz(:,1) = xyz(:,1)-xmean;
xyz(:,2) = xyz(:,2)-ymean;
xyz2 = [Cross(2).Snapline(:,1),Cross(2).Snapline(:,2),100*ones(size(Cross(2).Snapline,1),1)];
xmean1 = mean(xyz(:,1)); ymean1 = mean(xyz(:,2));
xyz2(:,1) = xyz2(:,1)-xmean1;
xyz2(:,2) = xyz2(:,2)-ymean1;
% calculate
psA = [xyz(:,1),xyz(:,2)];
psB = [xyz2(:,1),xyz2(:,2)];
sAsim = size(xyz,1); %point in A
sBsim = size(xyz2,1); %point in B
if sAsim < sBsim

    temp = psA;
    psA = psB;
    psB = temp;
end
sA = size(psA,1); %point in A
sB = size(psB,1); %point in B

matching = knnclassify(psA, psB, [1:1:sB]');
new_out = zeros(sA,2);
for step = 1:100
    for pti = 1:1:sA
        try
            if (new_out(pti,1)== psB(matching(pti),1)) && (new_out(pti,2) == psB(matching(pti),2))
                continue;   % decrease loop for corresponsedense point
                            % made lastest profile again
            end
            distX = pdist([psA(pti,1) ; psB(matching(pti),1)],'euclidean');
            distY = pdist([psA(pti,2) ; psB(matching(pti),2)],'euclidean');
            %% x and y is equal
            if (psA(pti,1) == psB(matching(pti),1)) && (psA(pti,2) == psB(matching(pti),2))
                new_out(pti,1) = psA(pti,1) ;
                new_out(pti,2) = psA(pti,2) ;
                
            %% x is different
            elseif psA(pti,1) ~= psB(matching(pti),1)
                if psA(pti,1) > psB(matching(pti),1)
                    nX =psA(pti,1)-(step*(distX/divide));
                
                        if nX < psB(matching(pti),1)
                        nX = psB(matching(pti),1);
                        end
                elseif psA(pti,1) < psB(matching(pti),1)
                    nX =psA(pti,1)+(step*(distX/divide));
                    if nX > psB(matching(pti),1)
                        nX = psB(matching(pti),1);
                    end
                end
                nY = interp1([psA(pti,1);psB(matching(pti),1)],[psA(pti,2);psB(matching(pti),2)],nX,'spline');
                new_out(pti,1) = nX;
                new_out(pti,2) = nY;
            %% x is not different
            else
                if psA(pti,2) > psB(matching(pti),2)
                    nX =psA(pti,2)-(step*(distY/divide));
                    if nX < psB(matching(pti),2)
                        nX = psB(matching(pti),2);
                    end
                elseif psA(pti,2) < psB(matching(pti),2)
                     nX =psA(pti,2)+(step*(distY/divide));
                    if nX > psB(matching(pti),2)
                        nX = psB(matching(pti),2);
                    end
                end
                nY = interp1([psA(pti,2);psB(matching(pti),2)],[psA(pti,1);psB(matching(pti),1)],nX,'spline');
                new_out(pti,1) = nY;
                new_out(pti,2) = nX;
            end
            
            if isnan(nY)
               disp(['ffff']) 
            end
        catch 
        end
        % =============== width condition ==================
        if pti == 1
            P1 = new_out(pti,1:2);
        elseif pti == ceil((sA-1)/2) % ======= experiment ======= 
            % half of shape 
            P2 = new_out(pti,1:2);
            D = pdist([P1;P2],'euclidean');
        end   
    end
    shape(step).XY = new_out;
end
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
[x_out,y_out,z_out] = rotationV3(B,BW,start,stop,Cross,shape);
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
