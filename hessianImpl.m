function [ BCell, CCell, E ] = hessianImpl(jacobCamCell, jacobPosCell, camInd, posInd, setting)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    jacobCamMat = zeros(setting(3)*2,setting(1)*9);
    jacobPosMat = zeros(setting(3)*2,setting(2)*3);
    for i = 1:setting(3)
        jacobCamMat(i*2-1:i*2,(camInd(i)-1)*9+1:camInd(i)*9) = jacobCamCell{i};
        jacobPosMat(i*2-1:i*2,(posInd(i)-1)*3+1:posInd(i)*3) = jacobPosCell{i};
    end
    E = jacobCamMat'*jacobPosMat;
    clear jacobCamMat jacobPosMat
    
    subJacCamCell = cell(1,setting(1));
    subJacPosCell = cell(1,setting(2));
    for i = 1:setting(3)
        subJacCamCell{camInd(i)}(end+1:end+2,:) = jacobCamCell{i};
        subJacPosCell{posInd(i)}(end+1:end+2,:) = jacobPosCell{i};
    end
    BCell = cellfun(@(x) x'*x, subJacCamCell, 'UniformOutput', false);
    CCell = cellfun(@(x) x'*x, subJacPosCell, 'UniformOutput', false);
end