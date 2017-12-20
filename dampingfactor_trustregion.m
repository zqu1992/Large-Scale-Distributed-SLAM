function [ mu_next, nu_next, delta2_next, update, f, stopLM ] = dampingfactor_trustregion(mu, nu, delta2, delta2_max, epsilon_3, epsilon_4, boundary, eta, stopLM, delta_m, delta_LM, f, argCam_next, argPos_next, pos2D_obs, camInd, posInd)
% Trust region radius (square radius: delta2) and damping factor updating in each step
%   Input:  mu: damping factor in this iteration step
%           mu_max: maximum damping factor
%           mu_min: minimum damping factor
%           delta2: square radius in this iteration step
%           delta2_max: maximum square radius (upper boundary)
%           boundary: if optimal increments step reaches delta2 (surface of ball) 
%           eta: threshold value for updating arguments
%           delta_m: increment value of model function
%           f: object function value in this iteration step
%           argCam_next: camera arguments for next iteration (with increments)
%           argPos_next: feature position arguments for next iteration (with increments)
%           pos2D_obs: observation of feature points (2D)
%           camInd: camera index of each observation 
%           posInd: feature position index of each observation
%   Output: mu_next: damping factor for next iteration step
%           delta2_next: square radius for next iteration step
%           update: if updating solution with optimal increments step in this iteration
%           f: object function value for next iteration step
    prjError_new = arrayfun(@projectionError, argCam_next(1,camInd),argCam_next(2,camInd),argCam_next(3,camInd),argCam_next(4,camInd),argCam_next(5,camInd),...
    argCam_next(6,camInd),argCam_next(7,camInd),argCam_next(8,camInd),argCam_next(9,camInd), argPos_next(1,posInd),argPos_next(2,posInd),argPos_next(3,posInd),...
    pos2D_obs(1,:),pos2D_obs(2,:), 'UniformOutput', false);
    prjError_new = reshape(cell2mat(prjError_new),1,[]);
    f_new = sum(bsxfun(@power, prjError_new, 2))/2;
    clear prjError prjError_new
    disp('f - f_new');
    disp(f - f_new);
    disp('delta_m');
    disp(delta_m);
    rho = (f - f_new)/delta_m;
    rho_LM = 2*(f - f_new)/delta_LM;
    % Choice 1:
    if (rho > eta)
        mu_next = max([mu/nu, 0.001]); % 4
        nu_next = 2;
    else
        mu_next = min([mu*nu, 10000]); % 2
        nu_next = 2;
    end
%     % Choice 2:
%     if (rho_LM > eta)
%         mu_next = mu * max(1/3, 1-(2*rho_LM-1)^3);
%         nu_next = 2;
%         %disp('epsilon_4');
%         %stopLM = stopLM || ((1 - epsilon_4) < sqrt(f_new/f));
%     else
%         mu_next = mu * nu; 
%         nu_next = nu * 2;
%     end
    if (rho<0.25)
        if (rho>0)
            delta2_next = delta2/16;
        else
            delta2_next = delta2;
        end
    else
        if (rho>0.75 && boundary)
            delta2_next = min([4*delta2, delta2_max]); 
        else
            delta2_next = delta2;
        end
    end
    update = false;
    if (f - f_new > 0)
        update = true;
        f = f_new;
    end
    %disp('epsilon_3');
    stopLM = stopLM || (f < epsilon_3);

end

