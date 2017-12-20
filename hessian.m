function [ B, C, E ] = hessian(jacobMat, setting)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    subJacCam = mat2cell(jacobMat(:,1:9*setting(1)), size(jacobMat,1), ones(1,setting(1))*9);
    diagBlockB = cellfun(@(x) x'*x, subJacCam, 'UniformOutput', false);
    B = blkdiag(diagBlockB{:});
    subJacPos = mat2cell(jacobMat(:,9*setting(1)+1:end), size(jacobMat,1), ones(1,setting(2))*3);
    diagBlockC = cellfun(@(x) x'*x, subJacPos, 'UniformOutput', false);
    C = blkdiag(diagBlockC{:});
    %B = origB + mu * diag(diag(origB));
    %C = origC + mu * diag(diag(origC));
    E = jacobMat(:,1:9*setting(1))'*jacobMat(:,9*setting(1)+1:end);
end

