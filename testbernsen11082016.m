% profile off;
function [] = testbernsen11082016(filename,pathname)
addpath('C:\Users\natar\Documents\MATLAB\masterProgram\A_master_Secondpaper\altmany-export_fig-8016f6a','-begin');
fig = figure;
%% bernsen 15
    T = 15;
    ws = [3 5 11 15 31]; 
    if iscell(filename)
        loop = size(filename,2);
    else
        loop = 1;
    end
for i = 1:1:loop
    %% manyfile or 1 file
    try
        tempname = filename{i};
    catch
        tempname = filename;
    end
    %% read path
    fich1 = fullfile(pathname,tempname);[~,name,~] = fileparts(tempname);
    
    %% if can't read skip it
    try
        rgb=imread(fich1);
    catch
        continue;
    end
    %% change to grayscale
    try
        im1 = rgb2gray(rgb);
    catch
        im1 = rgb;
    end
    
    for j = 1:1:size(ws,2)
        %% set path
        savingPath = strcat(pathname,'output\_result_TBernsenE2_T15\contour2\');
        Result_of = strcat('_result_BernsenE_T',num2str(T),'_ws',num2str(ws(j)));
        
        %% bernsen method
        [bwim1,~,~,localMask] = bernsenE2(im1,[ws(j) ws(j)],T);
        
        %% mask local contour
        BW2 = imfill(localMask,'holes');
        [B,L,N] = bwboundaries(BW2);
        set(fig,'name',strcat('compare',Result_of));
        
        %% image show
        imshow(BW2);        hold on;        warning('off', 'Images:initSize:adjustingMag');
        for m=1:length(B),
           boundary = B{m};
           plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
        end
        t(1) = text(2,10,strcat('number of contour',num2str(length(B))),'Color','blue');
        t(2) = text(2,50,strcat('window size ',num2str(ws(j)),'x',strcat(num2str(ws(j)))),'Color','blue');
        if ws(j) == 31
            % select largest contour
            [nrows,~] = cellfun(@size,B,'uni',false);
            [MVal,MInd] = max([nrows{:}]);
            boundary = B{MInd};
            plot(boundary(:,2), boundary(:,1), 'b','LineWidth',2);
            t(3) = text(2,30,num2str(MVal),'Color','blue');
        end
        
    %         t(3) =[];
    % imwrite(im1,strcat(savingPath,name,Result_of,'_An_origin','.jpg'));
    % imwrite(fig,strcat(savingPath,name,Result_of,'_detection','.jpg'));
%     print(fig,strcat(savingPath,name,Result_of,'_detection'),'-dpng');
%     export_fig strcat(savingPath,name,Result_of,'_detection') -jpg;
    imageData = export_fig;
    imwrite(imageData,strcat(savingPath,name,Result_of,'_detection','.jpg'));
    clf(fig);
    disp(strcat('loop Num: ',num2str(ws(j)),'_',tempname));
    t = [];
    if length(B) == 1
        break;
        
    end
   
    end
end