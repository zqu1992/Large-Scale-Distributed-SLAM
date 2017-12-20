function [ output_args ] = plotReprojection( argCam, argPos, camInd, posInd, pos2D_obs, setting)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    prjPos2D = cell2mat(arrayfun(@reprojectionCoord, argCam(1,camInd),argCam(2,camInd),argCam(3,camInd),argCam(4,camInd),argCam(5,camInd),...
        argCam(6,camInd),argCam(7,camInd),argCam(8,camInd),argCam(9,camInd), argPos(1,posInd),argPos(2,posInd),argPos(3,posInd),'UniformOutput', false));
    for i = 1:setting(1)
        figure(i)
        camInd_i = camInd == i;
        pos2D_obs_i = pos2D_obs(:, camInd_i);
        prjPos2D_i = prjPos2D(:, camInd_i);
        hold on
        plot(pos2D_obs_i(1,:), pos2D_obs_i(2,:), '*');
        plot(prjPos2D_i(1,:), prjPos2D_i(2,:), 'o');
        for j = 1:size(pos2D_obs_i,2)
            plot([pos2D_obs_i(1,j),prjPos2D_i(1,j)],[pos2D_obs_i(2,j),prjPos2D_i(2,j)],'b');
        end
        hold off
        
    end


end

