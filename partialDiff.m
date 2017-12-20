function [ dexey ] = partialDiff( phi1, phi2, phi3, trl1, trl2, trl3, f, k1, k2, pos3D1, pos3D2, pos3D3 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    phiMagSq = phi1^2 + phi2^2 + phi3^2;
    phiMag = sqrt(phiMagSq);
    cosPhiMag = cos(phiMag);
    sinPhiMag = sin(phiMag);

    dX2dphi1 = ((2*phi1*pos3D1 + phi2*pos3D2 + phi3*pos3D3)*(1-cosPhiMag) + phi1*cosPhiMag*(phi2*pos3D3 - phi3*pos3D2))/phiMagSq + ((2*phi1^3*pos3D1 + 2*phi1^2*phi2*pos3D2 + 2*phi1^2*phi3*pos3D3)*(cosPhiMag-1))/phiMagSq^2 + (phi1*sinPhiMag*(-phi2*pos3D3+phi3*pos3D2)+phi1^2*sinPhiMag*(phi2*pos3D2+phi3*pos3D3+phi1*pos3D1))/phiMagSq^(3/2) - (phi1*pos3D1*sinPhiMag)/phiMag;
    dX2dphi2 = ( pos3D3*sinPhiMag - phi2*pos3D1*sinPhiMag)/phiMag + (phi1*pos3D2*(1-cosPhiMag) + phi2*cosPhiMag*(phi2*pos3D3 - phi3*pos3D2))/phiMagSq + (2*phi1*phi2*(cosPhiMag-1)*(phi1*pos3D1 + phi2*pos3D2 + phi3*pos3D3))/phiMagSq^2 + (phi2*sinPhiMag*((phi1*phi3-phi2)*pos3D3 + (phi3+phi1*phi2)*pos3D2 + pos3D1*phi1^2))/phiMagSq^(3/2);
    dX2dphi3 = (-pos3D2*sinPhiMag - phi3*pos3D1*sinPhiMag)/phiMag + (phi1*pos3D3*(1-cosPhiMag) + phi3*cosPhiMag*(phi2*pos3D3 - phi3*pos3D2))/phiMagSq + (2*phi1*phi3*(cosPhiMag-1)*(phi1*pos3D1 + phi3*pos3D3 + phi2*pos3D2))/phiMagSq^2 + (phi3*sinPhiMag*(pos3D1*phi1^2 + (phi1*phi3-phi2)*pos3D3 + (phi1*phi2+phi3)*pos3D2))/phiMagSq^(3/2);
    dX2dpos3D1 = phi1^2/phiMagSq *(1-cosPhiMag) + cosPhiMag;
    dX2dpos3D2 = (phi1*phi2*(1-cosPhiMag))/phiMagSq - (phi3*sinPhiMag)/phiMag;
    dX2dpos3D3 = (phi1*phi3*(1-cosPhiMag))/phiMagSq + (phi2*sinPhiMag)/phiMag; 

    dY2dphi1 = (-pos3D3*sinPhiMag - phi1*pos3D2*sinPhiMag)/phiMag + (phi2*pos3D1*(1-cosPhiMag) + phi1*cosPhiMag*(-phi1*pos3D3 + phi3*pos3D1))/phiMagSq + (2*phi1*phi2*(cosPhiMag-1)*(phi2*pos3D2 + phi1*pos3D1 + phi3*pos3D3))/phiMagSq^2 + (phi1*sinPhiMag*((phi2*phi3+phi1)*pos3D3 + (-phi3+phi2*phi1)*pos3D1 + pos3D2*phi2^2))/phiMagSq^(3/2);
    dY2dphi2 = ((2*phi2*pos3D2 + phi1*pos3D1 + phi3*pos3D3)*(1-cosPhiMag) + phi2*cosPhiMag*(-phi1*pos3D3 + phi3*pos3D1))/phiMagSq + ((2*phi2^3*pos3D2 + 2*phi2^2*phi1*pos3D1 + 2*phi2^2*phi3*pos3D3)*(cosPhiMag-1))/phiMagSq^2 + (phi2*sinPhiMag*(phi1*pos3D3-phi3*pos3D1)+phi2^2*sinPhiMag*(phi1*pos3D1+phi2*pos3D2+phi3*pos3D3))/phiMagSq^(3/2) - (phi2*pos3D2*sinPhiMag)/phiMag;
    dY2dphi3 = ( pos3D1*sinPhiMag - phi3*pos3D2*sinPhiMag)/phiMag + (phi2*pos3D3*(1-cosPhiMag) + phi3*cosPhiMag*(-phi1*pos3D3 + phi3*pos3D1))/phiMagSq + (2*phi2*phi3*(cosPhiMag-1)*(phi1*pos3D1 + phi3*pos3D3 + phi2*pos3D2))/phiMagSq^2 + (phi3*sinPhiMag*(pos3D2*phi2^2 + (phi2*phi3+phi1)*pos3D3 + (phi1*phi2-phi3)*pos3D1))/phiMagSq^(3/2);
    dY2dpos3D1 = (phi2*phi1*(1-cosPhiMag))/phiMagSq + (phi3*sinPhiMag)/phiMag;
    dY2dpos3D2 = phi2^2/phiMagSq *(1-cosPhiMag) + cosPhiMag;
    dY2dpos3D3 = (phi2*phi3*(1-cosPhiMag))/phiMagSq - (phi1*sinPhiMag)/phiMag; 

    dZ2dphi1 = ( pos3D2*sinPhiMag - phi1*pos3D3*sinPhiMag)/phiMag + (phi3*pos3D1*(1-cosPhiMag) + phi1*cosPhiMag*(-phi2*pos3D1 + phi1*pos3D2))/phiMagSq + (2*phi3*phi1*(cosPhiMag-1)*(phi3*pos3D3 + phi1*pos3D1 + phi2*pos3D2))/phiMagSq^2 + (phi1*sinPhiMag*(pos3D3*phi3^2 + (phi3*phi1+phi2)*pos3D1 + (phi3*phi2-phi1)*pos3D2))/phiMagSq^(3/2);
    dZ2dphi2 = (-pos3D1*sinPhiMag - phi2*pos3D3*sinPhiMag)/phiMag + (phi3*pos3D2*(1-cosPhiMag) + phi2*cosPhiMag*(-phi2*pos3D1 + phi1*pos3D2))/phiMagSq + (2*phi3*phi2*(cosPhiMag-1)*(phi3*pos3D3 + phi2*pos3D2 + phi1*pos3D1))/phiMagSq^2 + (phi2*sinPhiMag*((phi3*phi1+phi2)*pos3D1 + (-phi1+phi3*phi2)*pos3D2 + pos3D3*phi3^2))/phiMagSq^(3/2);
    dZ2dphi3 = ((2*phi3*pos3D3 + phi2*pos3D2 + phi1*pos3D1)*(1-cosPhiMag) + phi3*cosPhiMag*(-phi2*pos3D1 + phi1*pos3D2))/phiMagSq + ((2*phi3^3*pos3D3 + 2*phi3^2*phi2*pos3D2 + 2*phi3^2*phi1*pos3D1)*(cosPhiMag-1))/phiMagSq^2 + (phi3*sinPhiMag*(phi2*pos3D1-phi1*pos3D2)+phi3^2*sinPhiMag*(phi2*pos3D2+phi1*pos3D1+phi3*pos3D3))/phiMagSq^(3/2) - (phi3*pos3D3*sinPhiMag)/phiMag;
    dZ2dpos3D1 = (phi3*phi1*(1-cosPhiMag))/phiMagSq - (phi2*sinPhiMag)/phiMag; 
    dZ2dpos3D2 = (phi3*phi2*(1-cosPhiMag))/phiMagSq + (phi1*sinPhiMag)/phiMag;
    dZ2dpos3D3 = phi3^2/phiMagSq *(1-cosPhiMag) + cosPhiMag;
    
    phi = [phi1;phi2;phi3];
    trl = [trl1;trl2;trl3];
    pos3D = [pos3D1;pos3D2;pos3D3];
    phiMat = phi*phi';
    pos3D_rot = (eye(3)-phiMat/phiMagSq)*cosPhiMag*pos3D + sinPhiMag/phiMag * cross(phi,pos3D)+phiMat/phiMagSq*pos3D;
    pos3D_cam = pos3D_rot + trl;
    p = -pos3D_cam(1:2)/pos3D_cam(3); % [px;py]
    pMagSq = p'*p;
    tmpDe = 1+k1*pMagSq+k2*pMagSq^2;
    dex2dpx = f + f*k1*(3*p(1)^2+p(2)^2)+f*k2*(5*p(1)^4+6*p(1)^2*p(2)^2+p(2)^4);
    dey2dpy = f + f*k1*(3*p(2)^2+p(1)^2)+f*k2*(5*p(2)^4+6*p(1)^2*p(2)^2+p(1)^4);
    dex2dpy = 2*f*p(1)*p(2)*(k1+2*k2*(p(1)^2+p(2)^2));
    dey2dpx = dex2dpy;   
    dpx2dX = -1/pos3D_cam(3);
    dpx2dZ = pos3D_cam(1)/pos3D_cam(3)^2;
    dpy2dY = dpx2dX;
    dpy2dZ = pos3D_cam(2)/pos3D_cam(3)^2;
    
    dex2df = p(1)*tmpDe;
    dey2df = p(2)*tmpDe;
    dex2dk1 = f*pMagSq*p(1);
    dey2dk1 = f*pMagSq*p(2);
    dex2dk2 = f*pMagSq^2*p(1);
    dey2dk2 = f*pMagSq^2*p(2);
    
    dex2dtrl1 = dex2dpx*dpx2dX;
    dex2dtrl2 = dex2dpy*dpy2dY;
    dex2dtrl3 = dex2dpx*dpx2dZ+dex2dpy*dpy2dZ;    
    dey2dtrl1 = dey2dpx*dpx2dX;
    dey2dtrl2 = dey2dpy*dpy2dY;
    dey2dtrl3 = dey2dpx*dpx2dZ+dey2dpy*dpy2dZ;
    
    dex2dphi1 = dex2dpx*(dpx2dX*dX2dphi1 + dpx2dZ*dZ2dphi1) + dex2dpy*(dpy2dY*dY2dphi1 + dpy2dZ*dZ2dphi1);
    dex2dphi2 = dex2dpx*(dpx2dX*dX2dphi2 + dpx2dZ*dZ2dphi2) + dex2dpy*(dpy2dY*dY2dphi2 + dpy2dZ*dZ2dphi2);
    dex2dphi3 = dex2dpx*(dpx2dX*dX2dphi3 + dpx2dZ*dZ2dphi3) + dex2dpy*(dpy2dY*dY2dphi3 + dpy2dZ*dZ2dphi3);
    dey2dphi1 = dey2dpx*(dpx2dX*dX2dphi1 + dpx2dZ*dZ2dphi1) + dey2dpy*(dpy2dY*dY2dphi1 + dpy2dZ*dZ2dphi1);
    dey2dphi2 = dey2dpx*(dpx2dX*dX2dphi2 + dpx2dZ*dZ2dphi2) + dey2dpy*(dpy2dY*dY2dphi2 + dpy2dZ*dZ2dphi2);
    dey2dphi3 = dey2dpx*(dpx2dX*dX2dphi3 + dpx2dZ*dZ2dphi3) + dey2dpy*(dpy2dY*dY2dphi3 + dpy2dZ*dZ2dphi3);
    
    dex2dpos3D1 = dex2dpx*(dpx2dX*dX2dpos3D1 + dpx2dZ*dZ2dpos3D1) + dex2dpy*(dpy2dY*dY2dpos3D1 + dpy2dZ*dZ2dpos3D1);
    dex2dpos3D2 = dex2dpx*(dpx2dX*dX2dpos3D2 + dpx2dZ*dZ2dpos3D2) + dex2dpy*(dpy2dY*dY2dpos3D2 + dpy2dZ*dZ2dpos3D2);
    dex2dpos3D3 = dex2dpx*(dpx2dX*dX2dpos3D3 + dpx2dZ*dZ2dpos3D3) + dex2dpy*(dpy2dY*dY2dpos3D3 + dpy2dZ*dZ2dpos3D3);
    dey2dpos3D1 = dey2dpx*(dpx2dX*dX2dpos3D1 + dpx2dZ*dZ2dpos3D1) + dey2dpy*(dpy2dY*dY2dpos3D1 + dpy2dZ*dZ2dpos3D1);
    dey2dpos3D2 = dey2dpx*(dpx2dX*dX2dpos3D2 + dpx2dZ*dZ2dpos3D2) + dey2dpy*(dpy2dY*dY2dpos3D2 + dpy2dZ*dZ2dpos3D2);
    dey2dpos3D3 = dey2dpx*(dpx2dX*dX2dpos3D3 + dpx2dZ*dZ2dpos3D3) + dey2dpy*(dpy2dY*dY2dpos3D3 + dpy2dZ*dZ2dpos3D3);

    dexey = [dex2dphi1,dex2dphi2,dex2dphi3,dex2dtrl1,dex2dtrl2,dex2dtrl3,dex2df,dex2dk1,dex2dk2,dex2dpos3D1,dex2dpos3D2,dex2dpos3D3;
             dey2dphi1,dey2dphi2,dey2dphi3,dey2dtrl1,dey2dtrl2,dey2dtrl3,dey2df,dey2dk1,dey2dk2,dey2dpos3D1,dey2dpos3D2,dey2dpos3D3];
end

