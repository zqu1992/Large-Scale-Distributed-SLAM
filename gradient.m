function [ grad ] = gradient( jacobMat, argCam, argPos, pos2D_obs, camInd, posInd)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    prjError = arrayfun(@projectionError, argCam(1,camInd),argCam(2,camInd),argCam(3,camInd),argCam(4,camInd),argCam(5,camInd),...
        argCam(6,camInd),argCam(7,camInd),argCam(8,camInd),argCam(9,camInd), argPos(1,posInd),argPos(2,posInd),argPos(3,posInd),...
        pos2D_obs(1,:),pos2D_obs(2,:), 'UniformOutput', false);
    prjError = reshape(cell2mat(prjError),[],1);
    grad = jacobMat'*prjError;
end

