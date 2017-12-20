
clear all
close all
clc

rng(1000);
%% Parameters Definition

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

%argCam = rand(9,setting(1));
%argPos = rand(3,setting(2));

argCam = argCamRes;
argPos = argPosRes;
%%
mode = 1;
numIter = round(setting(1)*9*0.7);
%numIter = 10;
delta2 = 25;
delta2_max = 100000;
eta = 0;
mu = 5;
mu_max = 100000;
mu_min = 0.01;
epsilon_1 = 1e-50; % stop LM requirement relevant to gradient norm
epsilon_2 = 1e-100; % stop LM requirement relavant to square increment step of arg
epsilon_3 = 1e-5; % stop LM requirement relavant to magnitude of f
epsilon_4 = 1e-500; % stop LM requirement relavant to reduction ratio of f     
numStep = 10000;
invSclVec = repmat(bsxfun(@power, 10, -[3,3,3,2,2,3,2,0,2]), 1, setting(1));
%% USE THIS
[arg_res, f_res, g_res, mu_res] = LM3(mu, delta2, delta2_max, epsilon_1, epsilon_2, epsilon_3, epsilon_4, eta, numStep, numIter, invSclVec, argCam, argPos, pos2D_obs, camInd, posInd, setting); 
% [arg, f_res, g_res, mu_res ] = gradientDesc(mu, numStep, argCam, argPos, pos2D_obs, camInd, posInd, setting);
