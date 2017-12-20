function [ arg, f_res, g_res, mu_res ] = gradientDesc(mu, numStep, argCam, argPos, pos2D_obs, camInd, posInd, setting )
% Levenberg-Marquardt algorithm with scaled trust-region and schur
% complement and PCG solver

    f_res = zeros(numStep,1);
    g_res = zeros(numStep,1);
    mu_res = zeros(numStep,1);
    k = 1;
    prjError = arrayfun(@projectionError, argCam(1,camInd),argCam(2,camInd),argCam(3,camInd),argCam(4,camInd),argCam(5,camInd),...
    argCam(6,camInd),argCam(7,camInd),argCam(8,camInd),argCam(9,camInd), argPos(1,posInd),argPos(2,posInd),argPos(3,posInd),...
    pos2D_obs(1,:),pos2D_obs(2,:), 'UniformOutput', false);
    prjError = reshape(cell2mat(prjError),1,[]);
    f = sum(bsxfun(@power, prjError, 2))/2;
    update = true;
    while(k<=numStep)
        if (update)
            [jacobCamCell, jacobPosCell] = jacobianImpl(argCam, argPos, setting, camInd, posInd);
            [ BCell, CCell, E ] = hessianImpl(jacobCamCell, jacobPosCell, camInd, posInd, setting);
            grad = gradientImpl(jacobCamCell, jacobPosCell, argCam, argPos, pos2D_obs, camInd, posInd, setting);         
            disp('norm g')
            disp(sum(bsxfun(@power,grad,2)))
        end
        argCam_next = argCam - mu*reshape(grad(1:setting(1)*9),9,[]);
        argPos_next = argPos - mu*reshape(grad(1+setting(1)*9:end),3,[]);
        [mu, update, f] = factorupdate_gradient(mu, f, argCam_next, argPos_next, pos2D_obs, camInd, posInd);
        
        if (update)
            argCam = argCam_next;
            argPos = argPos_next;
            f_res(k) = f;
            mu_res(k) = mu;
            g_res(k) = sum(abs(grad));
            k = k + 1;
        end
        disp('f')
        disp(f)
        
    end
    arg = [reshape(argCam,[],1);reshape(argPos,[],1)];
end

