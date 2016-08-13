%% TEST ALL CASE
clear all;close all;

[filename,pathname,index] = uigetfile({'*.*';'*.png';'*.jpg';'*.tif'}, 'Get Image file','MultiSelect','on');
fig = figure;
savingPath = strcat(pathname,'\output\_result_Tmean\');
Result_of = '_result_AdapTmean_k006_ws15';

ws = [3 5 11 15 31];
k = 0.06;
medianOrmean = 1;
    type = 'median';

for j = 1:size(ws,2)
    savingPath = strcat(pathname,'output\_result_T',type,'\');
    Result_of = strcat('_result_AdapT',type,'_k','006','_ws',num2str(ws(j)));
    

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
        rgb=imread(fich1);
        
        try
            im1 = rgb2gray(rgb);
        catch
            im1 = rgb;
        end
         bwim1=adaptivethreshold(im1,ws(j),k,medianOrmean);
         set(fig,'name',strcat('compare',Result_of));
         subplot(1,2,1);
         title('image');
         imshow(im1);
         subplot(1,2,2);
         title('adaptive threhold');
        imshow(bwim1);
% save file
savefig(fig,strcat(savingPath,name,Result_of,'.fig'));
imwrite(bwim1,strcat(savingPath,name,Result_of,'.jpeg'));
clf(fig);
% shutdown the warning

warning('off','MATLAB:print:SavingToDifferentName');
end
end
% online the warning
warning('on','MATLAB:print:SavingToDifferentName');

%% mean 
medianOrmean = 0;
type = 'mean';
for j = 1:size(ws,2)
    savingPath = strcat(pathname,'output\_result_T',type,'\');
    Result_of = strcat('_result_AdapT',type,'_k','006','_ws',num2str(ws(j)));
    

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
        rgb=imread(fich1);
        
        try
            im1 = rgb2gray(rgb);
        catch
            im1 = rgb;
        end
         bwim1=adaptivethreshold(im1,ws(j),k,medianOrmean);
         set(fig,'name',strcat('compare',Result_of));
         subplot(1,2,1);
         title('image');
         imshow(im1);
         subplot(1,2,2);
         title('adaptive threhold');
        imshow(bwim1);
% save file
savefig(fig,strcat(savingPath,name,Result_of,'.fig'));
imwrite(bwim1,strcat(savingPath,name,Result_of,'.jpeg'));
clf(fig);
% shutdown the warning

warning('off','MATLAB:print:SavingToDifferentName');
end
end
% online the warning
warning('on','MATLAB:print:SavingToDifferentName');

%% SAUVOLA
    k = 0.34;
    ws = [3 5 11 15 31];
for j = 1:size(ws,2)
    savingPath = strcat(pathname,'output\_result_TSauvola\');
    Result_of = strcat('_result_Sauvola_k034','_ws',num2str(ws(j)));
    
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
        rgb=imread(fich1);
        
        try
            im1 = rgb2gray(rgb);
        catch
            im1 = rgb;
        end
         bwim1 = sauvola(im1,[ws(j) ws(j)],k);
         set(fig,'name',strcat('compare',Result_of));
         subplot(1,2,1);
         title('image');
         imshow(im1);
         subplot(1,2,2);
         title('adaptive threhold');
        imshow(bwim1);
% save file
savefig(fig,strcat(savingPath,name,Result_of,'.fig'));
imwrite(bwim1,strcat(savingPath,name,Result_of,'.jpeg'));
clf(fig);
end
end

%% bernsen
    k = 128;
    ws = [3 5 11 15 31];
for j = 1:size(ws,2)
    savingPath = strcat(pathname,'output\_result_TBernsen\');
    Result_of = strcat('_result_Bernsen_T',num2str(k),'_ws',num2str(ws(j)));
    
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
        rgb=imread(fich1);
        
        try
            im1 = rgb2gray(rgb);
        catch
            im1 = rgb;
        end
         bwim1 = bernsen(im1,[ws(j) ws(j)],k);
         set(fig,'name',strcat('compare',Result_of));
         subplot(1,2,1);
         title('image');
         imshow(im1);
         subplot(1,2,2);
         title('adaptive threhold');
        imshow(bwim1);
% save file
savefig(fig,strcat(savingPath,name,Result_of,'.fig'));
imwrite(bwim1,strcat(savingPath,name,Result_of,'.jpeg'));
clf(fig);
end
end

%% 

%% new local T.R Singh’s 
    k = 0.06;
    ws = [3 5 11 15 31];
    
for j = 1:size(ws,2)
    savingPath = strcat(pathname,'output\_result_TAnewLocal\');
    Result_of = strcat('_result_AnewLocal_k006','_ws',num2str(ws(j)));
    
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
        rgb=imread(fich1);
        
        try
            im1 = rgb2gray(rgb);
        catch
            im1 = rgb;
        end
         bwim1 = NEWLOCAL(im1,[ws(j) ws(j)],k);
         set(fig,'name',strcat('compare',Result_of));
         subplot(1,2,1);
         title('image');
         imshow(im1);
         subplot(1,2,2);
         title('adaptive threhold');
        imshow(bwim1);
%% save file
savefig(fig,strcat(savingPath,name,Result_of,'.fig'));
imwrite(bwim1,strcat(savingPath,name,Result_of,'.jpeg'));
clf(fig);
end
end
