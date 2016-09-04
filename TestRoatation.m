% testing rotate
clear all;close all;
load TestDataHorizontal;
B = [B(:,2),B(:,1)];
[angle,pairMatching,midPoint] = calculateDegreeNPosition(B,start,stop);
for j = 1:100
profile(i) = CrossSection(1);
end
x = profile(1).shape(:,1);
y = profile(1).shape(:,2);
z = ones(1,size(x,2));
rotAxis = 2;
% start = position(1,:);
% stop = position(end,:);
for i = length(angle)
    profileScale(i) = norm(pairMatching(i,1:2) - pairMatching(i,3:4));
end
for i = 1:length(CrossSection)
    % rotation angle
    theta = angle(i)*pi/180;       % pi/3 radians = 60 degrees
    switch rotAxis
        case 1
            % 3D X
            R = [1          0           0; ...
                0           cos(theta)  -sin(theta); ...
                0           sin(theta)  cos(theta)];
        case 2
            % 3D Y
            R = [cos(theta) 0           sin(theta); ...
                0           1           0; ...
                -sin(theta) 0           cos(theta)];
        case 3
            % 3D Z
            R = [cos(theta) -sin(theta) 0; ...
                sin(theta)  cos(theta)  0; ...
                0           0           1];
    end
    % defined center of rotation
    v = [x;y;z];
%     s = regionprops(v,'centroid');
%     if start(1,1) < 0.5*size(I,1)
%         x_max = max((B(:,1)));
%     elseif start(1,1)>= 0.5*size(I,1)
%     end
    x_center = mean(x);
    y_center = mean(y);
    z_center = mean(z);
    center = repmat([x_center; y_center;z_center], 1, length(x));
%     Cscale = norm(pairMatching(i,1:2)-pairMatching(:,3:4));
%     Scale = abs(Cscale - profileScale);
%     if Scale == 0
%         Scale = 1 ;
%     end
    % rotate
    s = v - center;     % shift points in the plane so that the center of rotation is at the origin
    so = ((R)*Scale)*s;           % apply the rotation about the origin
    
    vo = so + center;   % shift again so the origin goes back to the desired center of rotation
    % this can be done in one line as:
    % vo = R*(v - center) + center
end










%% dump zone

% function [] = CrossRot(crosssection,angle,position)
% profile = crosssection(1);
% x = profile(1).XY(:,1);
% y = profile(1).XY(:,2);
% z = ones(1,size(x,2));
% start = position(1,:);
% stop = position(end,:);
% profileScale = norm(start(1,:) - start(2,:));
% for i = 1:length(crosssection)
%     % rotation angle
%     theta = angle(i)*pi/180;       % pi/3 radians = 60 degrees
%     switch rotAxis
%         case 1
%             % 3D X
%             R = [1          0           0; ...
%                 0           cos(theta)  -sin(theta); ...
%                 0           sin(theta)  cos(theta)];
%         case 2
%             % 3D Y
%             R = [cos(theta) 0           sin(theta); ...
%                 0           1           0; ...
%                 -sin(theta) 0           cos(theta)];
%         case 3
%             % 3D Z
%             R = [cos(theta) -sin(theta) 0; ...
%                 sin(theta)  cos(theta)  0; ...
%                 0           0           1];
%     end
%     % defined center
%     %     v = [x;y;z];
%     %     x_center = x(3);
%     %     y_center = y(3);
%     %     z_center = z(3);
%     v = position(i);
%     center = repmat([x_center; y_center;z_center], 1, length(x));
%     Cscale = norm(position(i)-position(i2));
%     Scale = abs(Cscale - profileScale);
%     if Scale == 0
%         Scale = 1 ;
%     end
%     % rotate
%     s = v - center;     % shift points in the plane so that the center of rotation is at the origin
%     so = ((R)*Scale)*s;           % apply the rotation about the origin
%     
%     vo = so + center;   % shift again so the origin goes back to the desired center of rotation
%     % this can be done in one line as:
%     % vo = R*(v - center) + center 
% end
% end