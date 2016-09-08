function [x_out,y_out,z_out] = rotationV3(B,BW,start,stop,CrossSection,shape)
%rotation by Arm & Benz

[angle,pairMatching,midPoint] = calculateDegreeNPosition(B,start,stop);

profile = CrossSection(1).Snapline;

x = profile(:,1)';
y = profile(:,2)';
z = zeros(size(x,2),1)';

temp = norm(pairMatching(1,1:2) - pairMatching(1,3:4));
for i = 1:length(angle)
    profileScale(i) = norm(pairMatching(i,1:2) - pairMatching(i,3:4));
    profileScale(i) = profileScale(i)./temp;
end

%% synth guide line
% define axis for start
x_out = cell(length(angle),1);
y_out = cell(length(angle),1);
z_out = cell(length(angle),1);
rotAxis = 2;

for i = 1:length(angle)
    % rotation angle
    theta = angle(i);
    theta = abs(theta-180);
    switch rotAxis
        case 1
            R = rotx(theta);
        case 2
            R = roty(theta);
        case 3
            R = rotz(theta);
    end
    
    % defined center of rotation
    v = [x;y;z];
    
    x_center = midPoint(i,1);
    y_center = 0;
    z_center = midPoint(i,2);
    
    center = repmat([x_center; y_center;z_center], 1, length(x));
    % rotate
    Scale = profileScale(i);
    v_p = v - repmat(mean(v, 2), 1, length(v));
    so = R*(v_p.*Scale);
    vo = so + center;   % shift again so the origin goes back to the desired center of rotation
    x_rotated = vo(1,:);
    y_rotated = vo(2,:);
    z_rotated = vo(3,:);
    hold on
    plot3(x_rotated, y_rotated, z_rotated, 'r-',x_center, y_center, z_center,'bo');
    x_out(i,:) = {x_rotated(:)};
    y_out(i,:) = {y_rotated(:)};
    z_out(i,:) = {z_rotated(:)};
    xlabel('x'), ylabel('y'), zlabel('z'), %text(x_center, y_center,z_center,num2str(angle(i)));
title('Perspective');
daspect([1 1 1]);
end
end