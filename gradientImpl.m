function [ grad ] = gradientImpl( jacobCamCell, jacobPosCell, argCam, argPos, pos2D_obs, camInd, posInd, setting)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    prjError = arrayfun(@projectionError, argCam(1,camInd),argCam(2,camInd),argCam(3,camInd),argCam(4,camInd),argCam(5,camInd),...
        argCam(6,camInd),argCam(7,camInd),argCam(8,camInd),argCam(9,camInd), argPos(1,posInd),argPos(2,posInd),argPos(3,posInd),...
        pos2D_obs(1,:),pos2D_obs(2,:), 'UniformOutput', false);
    
    subJacCamCell = cell(1,setting(1));
    subPrjCamCell = cell(1,setting(1));
    subJacPosCell = cell(1,setting(2));
    subPrjPosCell = cell(1,setting(2));
    for i = 1:setting(3)
        subJacCamCell{camInd(i)}(end+1:end+2,:) = jacobCamCell{i};
        subJacPosCell{posInd(i)}(end+1:end+2,:) = jacobPosCell{i};
        subPrjCamCell{camInd(i)}(end+1:end+2,:) = prjError{i};
        subPrjPosCell{posInd(i)}(end+1:end+2,:) = prjError{i};
    end
    grad = [reshape(cell2mat(cellfun(@(x,y) x'*y, subJacCamCell, subPrjCamCell, 'UniformOutput', false)),[],1);
        reshape(cell2mat(cellfun(@(x,y) x'*y, subJacPosCell, subPrjPosCell, 'UniformOutput', false)),[],1)];
end