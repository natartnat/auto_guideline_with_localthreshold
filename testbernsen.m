clear all;close all;
profile on;
% profile off;
[filename,pathname,~] = uigetfile({'*.*';'*.png';'*.jpg';'*.tif'}, 'Get Image file','MultiSelect','on');

fig = figure;
%% bernsen 15
    T = 15;
    a = 1;
    ws = [3 5 11 15 31];
for j = 1:size(ws,2)
    savingPath = strcat(pathname,'output\_result_TBernsenE2_T15\contour2\');
    Result_of = strcat('_result_BernsenE_T',num2str(T),'_ws',num2str(ws(j)));
    
    if iscell(filename)
        loop = size(filename,2);
    else
        loop = 1;
    end

for i = 1:1:loop
        try
            tempname = filename{i};
        catch
            tempname = filename;
        end
        fich1 = fullfile(pathname,tempname);[~,name,~] = fileparts(tempname);
        try
            rgb=imread(fich1);
        catch
            continue;
        end
        
        try
            im1 = rgb2gray(rgb);
        catch
            im1 = rgb;
        end
        set(fig,'name',strcat('compare',Result_of));
        [bwim1,mask,globalMask,localMask] = bernsenE2(im1,[ws(j) ws(j)],T);
        BW2 = imfill(localMask,'holes');
        [B,L,N] = bwboundaries(BW2);
        imshow(BW2);
        hold on;
        warning('off', 'Images:initSize:adjustingMag');
        for m=1:length(B),
           boundary = B{m};
           if(m > N)
             plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
           else
             plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
           end
        end
        t(1) = text(2,8,strcat('number of contour',num2str(length(B))));
        t(2) = text(2,20,strcat('window size ',num2str(ws(j)),'x',strcat(num2str(ws(j)))));
        if length(B) == 1
            filename{i} = [];
            
        elseif length(B) ~= 1 && ws(j) == 31
            % select largest contour
            [nrows,~] = cellfun(@size,B,'uni',false);
            [MVal,MInd] = max([nrows{:}]);
            boundary = B{MInd};
            plot(boundary(:,2), boundary(:,1), 'b','LineWidth',2);
            t(3) = text(2,30,num2str(MVal));
        end
        set(t(:),'color','g','fontw','bold','fonts',12);
%         t(3) =[];
% imwrite(im1,strcat(savingPath,name,Result_of,'_An_origin','.jpg'));
% imwrite(fig,strcat(savingPath,name,Result_of,'_detection','.jpg'));
print(fig,strcat(savingPath,name,Result_of,'_detection'),'-dpng');
clf(fig);
a = a+1;
disp(a);

end
end
profile viewer;
profile off;
