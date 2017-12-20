%%
clear all
clc
close all
load('optiRes2.mat')
%%
figure(1)
semilogy(1:100, f_res);
xlabel('step')
ylabel('objective function f')
figure(2)
plot(1:100, mu_res);
xlabel('step')
ylabel('damping factor \mu')
figure(3)
semilogy(1:100, g_res);
xlabel('step')
ylabel('norm of gradient')

%%
setting = load('setting.txt');
obs = load('observation.txt')';
arg = load('parameters.txt');
argCamRes = reshape(arg(1:setting(1)*9), 9,[]);
argPosRes = reshape(arg(1+setting(1)*9:end), 3,[]);
pos2D_obs = obs(3:4,:);
camInd = obs(1,:)+1;
posInd = obs(2,:)+1;
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
arg_res = [argCam(:);argPos(:)];
arg = [argCamRes(:);argPosRes(:)];
sum(abs(arg_res - arg))/sum(abs(arg))