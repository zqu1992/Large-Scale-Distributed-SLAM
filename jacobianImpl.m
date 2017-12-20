function [ jacobCamCell, jacobPosCell ] = jacobianImpl(argCam, argPos, setting, camInd, posInd)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    phi1 = argCam(1,camInd);
    phi2 = argCam(2,camInd);
    phi3 = argCam(3,camInd);
    trl1 = argCam(4,camInd);
    trl2 = argCam(5,camInd);
    trl3 = argCam(6,camInd);
    f = argCam(7,camInd);
    k1 = argCam(8,camInd);
    k2 = argCam(9,camInd);
    pos3D1 = argPos(1,posInd);
    pos3D2 = argPos(2,posInd);
    pos3D3 = argPos(3,posInd);
    jacobCell = arrayfun(@partialDiff, phi1, phi2, phi3, trl1, trl2, trl3, f, k1, k2, pos3D1, pos3D2, pos3D3, 'UniformOutput', false);
    jacobCamCell = cell(1,setting(3));
    jacobPosCell = cell(1,setting(3));
    for i = 1:setting(3)
        jacobCamCell{i} = jacobCell{i}(:,1:9);
        jacobPosCell{i} = jacobCell{i}(:,10:12);
    end

end
