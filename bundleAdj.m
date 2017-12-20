clear all
close all
clc
setting = load('setting.txt');
obs = load('observation.txt')';
arg = load('parameters.txt');
%% Load Data
argCamRes = reshape(arg(1:setting(1)*9), 9,[]);
argPosRes = reshape(arg(1+setting(1)*9:end), 3,[]);
pos2D_obs = obs(3:4,:);
camInd = obs(1,:)+1;
posInd = obs(2,:)+1;
%% Reduce Dataset Scale
setting(1) = 10;
tmpCamInd = camInd > setting(1);
argCamRes(:,setting(1)+1:end) = [];
camInd(tmpCamInd) = [];
posInd(tmpCamInd) = [];
pos2D_obs(:,tmpCamInd) = [];
tbl = tabulate(posInd);
validPosInd = sort(tbl(tbl(:,2)>1,1));
tmpPosInd = ismember(posInd,validPosInd);
posInd = posInd(tmpPosInd);
camInd = camInd(tmpPosInd);
pos2D_obs = pos2D_obs(:,tmpPosInd);
argPosRes = argPosRes(:, validPosInd);
setting(2) = length(validPosInd);
for i = 1:length(validPosInd)
    posInd(posInd == validPosInd(i)) = i;   
end
setting(3) = length(camInd);
%% Reprojection Position Demonstration
% 
load('optiRes_choice1.mat');
argCamRes = reshape(arg_res(1:setting(1)*9,4),9,[]);
argPosRes = reshape(arg_res(1+setting(1)*9:end,4),3,[]);
plotReprojection(argCamRes, argPosRes, camInd, posInd, pos2D_obs, setting);

% plotReprojection(argCamRes, argPosRes, camInd, posInd, pos2D_obs, setting);
%% Reprojection Error
% calculate the reprojection error based on the optimal solution from dataset 
prjError = arrayfun(@projectionError, argCamRes(1,camInd),argCamRes(2,camInd),argCamRes(3,camInd),argCamRes(4,camInd),argCamRes(5,camInd),...
    argCamRes(6,camInd),argCamRes(7,camInd),argCamRes(8,camInd),argCamRes(9,camInd), argPosRes(1,posInd),argPosRes(2,posInd),argPosRes(3,posInd),...
    pos2D_obs(1,:),pos2D_obs(2,:), 'UniformOutput', false);
prjError = reshape(cell2mat(prjError),1,[]);
histogram(prjError)
title('reprojection error F(x)')
sum(bsxfun(@power, prjError, 2))/2
sum(abs(prjError))/setting(3)
%% Analytical Jacobian Matrix
jacobMat = jacobian(argCamRes, argPosRes, setting, camInd, posInd);
%% Numerical Jacobian Matrix
delta = [1e-8,1e-8,1e-8,1e-5,1e-5,1e-8,1e0,1e-3,1e-3,1e-7,1e-7,1e-8]; % 0.7204 2.9383 
%delta = delta / 2;
%delta = [1e-10,1e-9,1e-7,1e-5,1e-5,1e-8,1e0,1e-2,1e-2,1e-7,1e-7,1e-8]; % 8.1741 2.9364
jacobMatNumerical = jacobianNumerical(argCamRes, argPosRes, pos2D_obs, setting, camInd, posInd, delta);
sum(sum(abs(bsxfun(@minus, jacobMat, jacobMatNumerical))))/sum(sum(abs(jacobMat)))
%norm(bsxfun(@minus, jacobMat, jacobMatNumerical))
%% Verification
tau = -5:30;
tau = bsxfun(@power, 2, -tau);
deltaCam = rand(size(argCamRes));
deltaPos = rand(size(argPosRes));
smallTerm = zeros(size(tau));
jacobDelta = jacobMat * [deltaCam(:);deltaPos(:)];
for i = 1:length(tau)
    argCamResNew = argCamRes + tau(i)*deltaCam;
    argPosResNew = argPosRes + tau(i)*deltaPos;
    prjErrorNew = arrayfun(@projectionError, argCamResNew(1,camInd),argCamResNew(2,camInd),argCamResNew(3,camInd),argCamResNew(4,camInd),argCamResNew(5,camInd),...
    argCamResNew(6,camInd),argCamResNew(7,camInd),argCamResNew(8,camInd),argCamResNew(9,camInd), argPosResNew(1,posInd),argPosResNew(2,posInd),argPosResNew(3,posInd),...
    pos2D_obs(1,:),pos2D_obs(2,:), 'UniformOutput', false);
    prjErrorNew = reshape(cell2mat(prjErrorNew),1,[]);
    smallTerm(i) = sqrt(sum(bsxfun(@power, (prjErrorNew-prjError)' - tau(i)*jacobDelta, 2))); 
end
plot(1:length(tau),log(1./smallTerm)/log(2));
grid on
%% Optimal Numerical Step Determination 
delta = [1e-8,1e-8,1e-8,1e-5,1e-5,1e-8,1e0,1e-3,1e-3,1e-7,1e-7,1e-8];
minNorm = 1e10;
for i = 1:12
    delta_i = delta(i);
    for j = 0:-15
        delta(i) = 10^(-j);
        jacobMatNumerical = jacobianNumerical(argCamRes, argPosRes, pos2D_obs, setting, camInd, posInd, delta);
        norm_now = sum(sum(abs(bsxfun(@minus, jacobMat, jacobMatNumerical))));
        if(norm_now < minNorm)
            minNorm = norm_now;
            delta_i = delta(i);
        end
    end
    delta(i) = delta_i;
end
%% Hessian Matrix
[B, C, E] = hessian(jacobMat, setting);
hessianMat1 = [B,E;E',C];
%det(hessianMat1)
%% Implicity Jacobian Matrix
[jacobCamCell, jacobPosCell] = jacobianImpl(argCamRes, argPosRes, setting, camInd, posInd);
%% Implicity Hessian Matrix
[BCell, CCell, E] = hessianImpl(jacobCamCell, jacobPosCell, camInd, posInd, setting);
B = blkdiag(BCell{:});
C = blkdiag(CCell{:});
hessianMat2 = [B,E;E',C];
sum(sum(abs(bsxfun(@minus,hessianMat1, hessianMat2))))
%% Verify Correctness of Hessian Matrix
hessianMat3 = jacobMat'*jacobMat;
sum(sum(abs(bsxfun(@minus,hessianMat1, hessianMat3))))
%% Gradient Array
grad1 = gradient(jacobMat, argCamRes, argPosRes, pos2D_obs, camInd, posInd);
%% Implicity Gradient Array
grad2 = gradientImpl(jacobCamCell, jacobPosCell, argCamRes, argPosRes, pos2D_obs, camInd, posInd, setting);
sum(bsxfun(@power, grad2, 2))
histogram(abs(grad2))
title('grad norm')
sum(sum(abs(bsxfun(@minus,grad1, grad2))))
