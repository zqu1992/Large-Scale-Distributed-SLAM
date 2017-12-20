function [ mu_next, update, f] = factorupdate_gradient(mu, f, argCam_next, argPos_next, pos2D_obs, camInd, posInd)
    prjError_new = arrayfun(@projectionError, argCam_next(1,camInd),argCam_next(2,camInd),argCam_next(3,camInd),argCam_next(4,camInd),argCam_next(5,camInd),...
    argCam_next(6,camInd),argCam_next(7,camInd),argCam_next(8,camInd),argCam_next(9,camInd), argPos_next(1,posInd),argPos_next(2,posInd),argPos_next(3,posInd),...
    pos2D_obs(1,:),pos2D_obs(2,:), 'UniformOutput', false);
    prjError_new = reshape(cell2mat(prjError_new),1,[]);
    f_new = sum(bsxfun(@power, prjError_new, 2))/2;
    clear prjError prjError_new
    rho = f - f_new;
    disp('f-f_new')
    disp(rho)
    if (rho > 0)
        mu_next = mu*2; % 4
    else
        mu_next = mu/2; % 2
    end
    update = false;
    if (rho > 0)
        update = true;
        f = f_new;
    end
end

