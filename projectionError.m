function [ error ] = projectionError( phi1, phi2, phi3, trl1, trl2, trl3, f, k1, k2, pos3D1, pos3D2, pos3D3, pos2D_obs1, pos2D_obs2)
% Single reprojection error corresponding to one single observation 
%   Input:  phi1, phi2, phi3: camera pose parameters (three rotation angles)
%           trl1, trl2, trl3: camera pose parameters (three translation coordinates) 
%           f: camera pose parameters (focal length)
%           k1, k2: camera pose parameters (two distortion calibration parameters)
%           pos3D1, pos3D2, pos3D3: feature position coordinates (3D)
%           pos2D_obs1, pos2D_obs2: observation point coordinates (2D)
%   Output: error: single reprojection error (2D)  

    phi = [phi1;phi2;phi3];
    trl = [trl1;trl2;trl3];
    pos3D = [pos3D1;pos3D2;pos3D3];
    pos2D_obs = [pos2D_obs1;pos2D_obs2];
    phiMat = phi*phi';
    phiDot = phi'*phi;
    % @demmeln FIXME: if phi is small, use Taylor expansion @Zhongnan Yes,
    % I Will try it.
    phiAbs = sqrt(phiDot);
    pos3D_rot = (eye(3)-phiMat/phiDot)*cos(phiAbs)*pos3D + sin(phiAbs)/phiAbs* cross(phi,pos3D)+phiMat/phiDot*pos3D;
    pos3D_trl = pos3D_rot + trl;
    % FIXME @demmeln: what z=0? @Zhongnan Yes, actually in a frame, there are only two coordindates. 
    pos2D_div = -pos3D_trl(1:2)/pos3D_trl(3);
    pos2DSq = pos2D_div'*pos2D_div;
    pos2D = f * (1+k1*pos2DSq+k2*pos2DSq^2) * pos2D_div;
    error = pos2D-pos2D_obs;
    
end

