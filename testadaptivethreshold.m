clear all;close all;

[filename,pathname,index] = uigetfile({'*.*';'*.png';'*.jpg';'*.tif'}, 'Get Image file','MultiSelect','on');

% savingPath = strcat(pathname,'\output\_result_Tmean\');
% Result_of = '_result_AdapTmean_k006_ws15';
savingPath = strcat(pathname,'\output\_result_Tmedian\');
Result_of = '_result_AdapTmedian_k006_ws15';
ws = 15;
k = 0.06;
medianOrmean = 1;

fig = figure;
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
         bwim1=adaptivethreshold(im1,ws,k,medianOrmean);
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
%% shutdown the warning

warning('off','MATLAB:print:SavingToDifferentName');
end
clear tempname;
close all;
%% online the warning
warning('on','MATLAB:print:SavingToDifferentName');
 % im1=imread('page.png');
% bwim1=adaptivethreshold(im1,11,0.03,0);
% subplot(2,2,1);
% imshow(im1);
% subplot(2,2,2);
% imshow(bwim1);

% im2=imread('tshape.png');
% bwim2=adaptivethreshold(im2,15,0.02,0);
% subplot(2,2,3);
% imshow(im2);
% subplot(2,2,4);
% imshow(bwim2);