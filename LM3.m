function [ arg, f_res, g_res, mu_res ] = LM3(mu, delta2, delta2_max, epsilon_1, epsilon_2, epsilon_3, epsilon_4, eta, numStep, numIter, invSclVec, argCam, argPos, pos2D_obs, camInd, posInd, setting )
% Levenberg-Marquardt algorithm with scaled trust-region and schur
% complement and PCG solver
    nu = 2;
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
    stopLM = false;
    while(k<=numStep && ~stopLM )
        if (update)
            [jacobCamCell, jacobPosCell] = jacobianImpl(argCam, argPos, setting, camInd, posInd);
            [ BCell, CCell, E ] = hessianImpl(jacobCamCell, jacobPosCell, camInd, posInd, setting);
            grad = gradientImpl(jacobCamCell, jacobPosCell, argCam, argPos, pos2D_obs, camInd, posInd, setting);
        end
        if (max(abs(grad))<epsilon_1)
            disp('epsilon_1')
            break;
        end
        % Disable trust region boundary
        % delta2 = 1e100;
        [ P, delta_m, boundary ] = PCG_SC_TR_SCALED(mu, delta2, numIter, invSclVec, BCell, CCell, E, grad(1:setting(1)*9), grad(1+setting(1)*9:end), setting);
        mag2_P = sum(bsxfun(@power, P, 2));
        if (mean(bsxfun(@rdivide, bsxfun(@power, P, 2), bsxfun(@power, [argCam(:);argPos(:)], 2)))<epsilon_2)
        %if (mag2_P <= epsilon_2 * sum(bsxfun(@power, [argCam(:);argPos(:)], 2))) % differ with the condition in the paper
            disp('epsilon_2')
            break;
        else
            argCam_next = argCam + reshape(P(1:setting(1)*9),9,[]);
            argPos_next = argPos + reshape(P(1+setting(1)*9:end),3,[]);  
            delta_LM = P'*-grad + mag2_P * mu;
            [ mu, nu, delta2, update, f, stopLM ] = dampingfactor_trustregion(mu, nu, delta2, delta2_max, epsilon_3, epsilon_4, boundary, eta, stopLM, delta_m, delta_LM, f, argCam_next, argPos_next, pos2D_obs, camInd, posInd);
            %disp(delta_m);
            if (update)
                argCam = argCam_next;
                argPos = argPos_next;
                f_res(k) = f;
                mu_res(k) = mu;
                g_res(k) = sum(abs(grad));
                k = k + 1;
            end
        end
        disp('f')
        disp(f)
        
    end
    arg = [reshape(argCam,[],1);reshape(argPos,[],1)];
end

