function [ jacobMat ] = jacobianNumerical( argCam, argPos, pos2Dgt, setting, camInd, posInd, delta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    jacobMat = zeros(setting(3)*2,setting(1)*9+setting(2)*3);
    
    prjError0 = cell2mat(arrayfun(@projectionError, argCam(1,camInd),argCam(2,camInd),argCam(3,camInd),argCam(4,camInd),argCam(5,camInd),...
        argCam(6,camInd),argCam(7,camInd),argCam(8,camInd),argCam(9,camInd), argPos(1,posInd),argPos(2,posInd),argPos(3,posInd),...
        pos2Dgt(1,:),pos2Dgt(2,:), 'UniformOutput', false));
    for i = 1:setting(3)
        argDelta = bsxfun(@plus, [argCam(:,camInd(i));argPos(:,posInd(i))], diag(delta));
        pos2DgtTmp = repmat(pos2Dgt(:,i),1,12);
        prjErrorDelta = cell2mat(arrayfun(@projectionError, argDelta(1,:),argDelta(2,:),argDelta(3,:),argDelta(4,:),...
            argDelta(5,:),argDelta(6,:),argDelta(7,:),argDelta(8,:),argDelta(9,:),argDelta(10,:),argDelta(11,:),...
            argDelta(12,:),pos2DgtTmp(1,:),pos2DgtTmp(2,:), 'UniformOutput', false));
        grad = bsxfun(@rdivide, bsxfun(@minus, prjErrorDelta, prjError0(:,i)), delta);
        jacobMat(2*i-1:2*i,(camInd(i)-1)*9+1:camInd(i)*9) = grad(:,1:9);
        jacobMat(2*i-1:2*i,setting(1)*9+(posInd(i)-1)*3+1:setting(1)*9+posInd(i)*3) = grad(:,10:12);
    end

end

