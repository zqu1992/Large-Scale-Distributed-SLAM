function [ P, delta_m, boundary ] = PCG_SC_TR_SCALED(mu, delta2, numIter, invSclVec, BCell, CCell, E, gradCam, gradPos, setting)
% Preconditioned conjugate gradient based on scaled trust region radius,
% and Schur Complement
%   Input:  BCell: Matrix B in cell type
%           CCell: Matrix B in cell type
%           E: Matrix E in matrix type
%           gradCam: gradient array corresponding to camera 
%           gradPos: gradient array corresponding to feature position
%           numIter: maximun iteration number of PCG solver
%           mu: damping factor of Levenberg-Marquardt
%           delta2: square radius of trust region (ball)
%           invSclVec: inverse scaling matrix (only diagonal element)
%           setting: dataset characteristics
%   Output: P: optimal increment step of whole arg space in trust region
%           delta_m: increment value of model function
%           boundary: if P is located on the boundary of trust region (surface)

%%%%%%% Deficiency of provement: if delta2 for cam args can be updated
%%%%%%% based on rho which is calculate based on whole args space

%%%%%%% Deficiency of provement: if delta2 and mu (damping factor) can be
%%%%%%% both updated based on rho at the same time, or if this updating
%%%%%%% process suitable


% Levenberg-Marquardt damping factor: H = J'*J + mu * D'*D = [B, E; E' C]
% H * [Py; Pz] = -[v ; w]
% Choice 1: D'*D = diag(J'*J)
% CCell = cellfun(@(x) x+mu*diag(diag(x)), CCell, 'UniformOutput', false);
% BCell = cellfun(@(x) x+mu*diag(diag(x)), BCell, 'UniformOutput', false);
% Choice 2: D'*D = I
CCell = cellfun(@(x) x+mu*eye(size(x)), CCell, 'UniformOutput', false);
BCell = cellfun(@(x) x+mu*eye(size(x)), BCell, 'UniformOutput', false);
% Choice 3: no damping factor (Gaussian-Newton)
% CCell = CCell;
% BCell = BCell;

% Schur Complement: S*Py = -u  
% with S = (B - E*inv(C)*E'), u = v - E*inv(C)*w
CCell = cellfun(@condAdjust, CCell, 'UniformOutput', false);
invC = cellfun(@(x) inv(x), CCell, 'UniformOutput', false);
invC = sparse(blkdiag(invC{:}));
B = sparse(blkdiag(BCell{:}));
u = invC * gradPos;
rem = gradPos' * u * 0.5;
u = gradCam - E * u;
if (isnan(rem))
    disp('wrong');
end
% Scaling: u_scl = inv(sclMat)*u, S_scl = inv(sclMat) * S * inv(sclMat),
% Py_scl = sclMat * Py
u = diag(invSclVec) * u;

mag2_u = sqrt(sum(bsxfun(@power, u, 2)));
epsilon2 = min([0.25, sqrt(mag2_u)])*mag2_u; % superlinear

% Preconditioner: M 
% Choice 1: M = B
% invMCell = cellfun(@(x) inv(x), BCell, 'UniformOutput', false);
% Choice 2: M = blkdiag(S)
% MCell = reshape(cellfun(@(x,y) x*y', mat2cell(E*invC, zeros(1, setting(1))+9),...
%     mat2cell(E, zeros(1, setting(1))+9), 'UniformOutput', false),1,setting(1));
% MCell = cellfun(@(x,y) x-y, BCell, MCell, 'UniformOutput', false);
% invMCell = cellfun(@(x) inv(x), MCell, 'UniformOutput', false);
% Choice 3: no Preconditioner, M = I
invMCell = cell(1,setting(1));
for i = 1:setting(1)
    invMCell{i}= eye(9);
end
% Judge if inv(M) has full rank
% if (rank(blkdiag(invMCell{:}))<90)
%     disp('Wrong')
% end

% Initial value of PCG 
p = zeros(setting(1)*9,1);
r = u;
q = blkdiag(invMCell{:})*r;
d = -q;
boundary = false;

% Preconditioned Conjugate Gradient: solve optimal increment step Py with preconditioner M
% PCG maximum excuted numIter times 
for j = 1:numIter
    scld = diag(invSclVec) * d;
    Sd = E'*scld;
    Sd = invC * Sd;
    Sd = B*scld - E*Sd;
    Sd = diag(invSclVec) * Sd;
    alpha = r'*q/(d'*Sd);
    p = p + alpha * d;
    mag2_p = sum(bsxfun(@power, p, 2));
    % Judge if solution exceeds ball radius of trust region
    if (mag2_p >= delta2)
        dot_pd = 2 * sum(bsxfun(@times, p, d));
        mag2_d = sum(bsxfun(@power, d, 2));
        tau = (-dot_pd + sqrt(dot_pd^2 - 4*mag2_d*(mag2_p - delta2)))/(2*mag2_d);
        Py = p + tau * d;
        sclPy = diag(invSclVec) * Py;
        Sp = E'*sclPy;
        Pz = -gradPos - Sp;
        Sp = invC* Sp;
        Sp = B*sclPy - E*Sp;
        Sp = diag(invSclVec) * Sp;
        % Find increment value of model function
        delta_m = -u'*Py - Py'*Sp*0.5 + rem;
        % Update whole increment step
        P = [sclPy; invC*Pz];
        boundary = true;
        return
    end
    r_next = r + alpha * Sd;
    % Judge if PCG has find the optimal increment step, and jump out of PCG
    % loop
    if (sum(bsxfun(@power, r_next, 2)) < epsilon2)
        Py = p;
        sclPy = diag(invSclVec) * Py;
        Sp = E'*sclPy;
        Pz = -gradPos - Sp;
        Sp = invC * Sp;
        Sp = B*sclPy - E*Sp;
        Sp = diag(invSclVec) * Sp;
        % Find increment value of model function
        delta_m = -u'*Py - Py'*Sp*0.5 + rem;
        % Update whole increment step
        P = [sclPy; invC*Pz];
        return
    end
    q_next = blkdiag(invMCell{:})*r_next;
    beta = r_next'*q_next/(r'*q);
    d = -q_next + beta * d;
    q = q_next;
    r = r_next;
end
scld = diag(invSclVec) * d;
Sd = E'*scld;
Sd = invC * Sd;
Sd = B*scld - E*Sd;
Sd = diag(invSclVec) * Sd;
alpha = r'*q/(d'*Sd);
p = p + alpha * d;
Py = p;
sclPy = diag(invSclVec) * Py;
Sp = E'*sclPy;
Pz = -gradPos - Sp;
Sp = invC * Sp;
Sp = B*sclPy - E*Sp;
Sp = diag(invSclVec) * Sp;
% Find increment value of model function
delta_m = -u'*Py - Py'*Sp*0.5 + rem;
% Update whole increment step
P = [sclPy; invC*Pz];
return

end

