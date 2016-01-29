
% [u, errAll] = [u,errAll] = PICCS_fMRI(RAll,fAll, dimIm,uref,muu, lambda,
%                                   alpha, gamma, nBreg, imTrueAll, h4)
%
% Inputs:
%
% RAll      = undersampling matrix, same size as fAll
% fAll      = 2D+time data, which corresponds to fft2(imTrueAll)+noise [complex]
% dimIm     = Nx*Ny*frames
% uref      = prior image [complex]
% muu       = parameter weighting the data fidelity term, use mu=1
% lambda    = parameter weighting the constraints, use lambda=1
% alpha     = parameter weighting the prior term, use=0.95 and tune depending of the problem
% gamma     = parameter to improve the conditioning, use between mu/100 and
% nBreg     = number of (outer) iterations
% imTrueAll = target image to compute the error at each iteration [abs]
% h4        = figure handle for display
%
% Outputs:
%
% u         = reconstructed image, size dimIm
% errAll    = Relative solution error norm at each iteration for all frames
%
% -------------------------------------------------------------------------
% PICCS is solved using a modification of Goldstein'n code mrics.m
% downloaded from
% (http://www.ece.rice.edu/~tag7/Tom_Goldstein/Split_Bregman.html), see Tom
% Goldstein and Stanley Osher. The Split Bregman Method for L1-Regularized
% Problems. SIAM J. Imaging Sci., 2(2), 323–343.
%
% -------------------------------------------------------------------------
% If you use this code, please, cite the following paper: Chavarrias, C.,
% Abascal, J.F., Montesinos, P. and Desco, M.,"Exploitation of temporal
% redundancy in compressed sensing reconstruction of fMRI studies with a
% prior-based algorithm (PICCS)". Med Phys,  42(7): p. 3814 (2015). DOI:
% http://dx.doi.org/10.1118/1.4921365
%
% Code downloaded from repository:
% https://github.com/HGGM-LIM/Split_Bregman_PICCS_fMRI
%
%
%
%
% -------------------------------------------------------------------------
% Cristina Chavarrías, Juan FPJ Abascal
% Departamento de Bioingeniería e Ingeniería Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% cristina.chavarrias@gmail.com, juanabascal78@gmail.com, desco@hggm.es


function [u,errAll] = PICCS_fMRI(RAll,fAll, dimIm,uref,muu, lambda,alpha, gamma, nBreg, imTrueAll, h4)


rows        =   dimIm(1);
cols        =   dimIm(2);
numTime     =   dimIm(3);

errAll      =   zeros(nBreg,numTime);


% normalize the data so that standard parameter values work
normFactor  =   getNormalizationFactor(RAll(:,:,1),fAll(:,:,1));
fAll        =   normFactor*fAll;

uref        =   uref*normFactor;
imTrueAll   =   imTrueAll*normFactor;


% Reserve memory for the auxillary variables
f0All       =   fAll;
u           =   zeros(dimIm);
x           =   zeros(dimIm);
y           =   zeros(dimIm);
bx          =   zeros(dimIm);
by          =   zeros(dimIm);

xe          =   zeros(dimIm);
ye          =   zeros(dimIm);

cx          =   zeros(dimIm);
cy          =   zeros(dimIm);

dxref       =   Dx(uref);
dyref       =   Dy(uref);

if size(uref,3) == 1
    dxref    =  repmat(dxref,[1 1 numTime]);
    dyref    =  repmat(dyref,[1 1 numTime]);
end
uBest       =   u;
errBest     =   inf;



% Build Kernels
% Spatial Hessian
uker        =   zeros(rows,cols);
uker(1,1)   =   4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
uker        =   (lambda + lambda)*fft2(uker);
uker        =   repmat(uker,[1 1 numTime]);
uker        =   uker + gamma + muu*RAll;

murf        =   ifft2(muu*fAll);


h   = waitbar(0);
%  Do the reconstruction
for outer = 1:nBreg;
    
    % update u
    rhs    =   murf+lambda*(Dxt(x-bx)+Dyt(y-by))+ ...
        + lambda*(Dxt(xe+dxref-cx)+Dyt(ye+dyref-cy))+gamma*u;
    u      =   ifft2(fft2(rhs)./uker);
    
    % update x and y
    dx      =   Dx(u);
    dy      =   Dy(u);
    [x,y]   =   shrink2(dx+bx, dy+by, (1-alpha)/lambda);
    [xe,ye] =   shrink2(dx-dxref+cx, dy-dyref+cy,alpha/lambda);
    
    
    % update bregman parameters
    bx      =   bx+dx-x;
    by      =   by+dy-y;
    
    cx      =   cx+dx-dxref-xe;
    cy      =   cy+dy-dyref-ye;
    
    fForw   =   RAll.*fft2(u);
    fAll    =   fAll + f0All - fForw;
    murf    =   ifft2(muu*fAll);
    
    
    % Compute the error
    uNorm           =   reshape(abs(imTrueAll),rows*cols,numTime);
    uNorm           =   sqrt(sum(uNorm.*uNorm));
    errThis         =   reshape(abs(imTrueAll-u),rows*cols,numTime);
    errThis         =   sqrt(sum(errThis.*errThis))./uNorm;
    errAll(outer,:) =   errThis;
    
    
    if (mean(errThis) <= errBest)
        uBest       = u;
        errBest     = mean(errThis);
    end %errThis
    
    % Displays the reconstructed first volume of the series and plots 
    % 1)the average error norm 
    % 2)the time course of the pixel [22,20] in the reconstructed fMRI series
    if any([outer <6 outer==10 mod(outer,50)==0])
        figure(h4),
        subplot(2,2,1); imagesc(abs(u(:,:,1))), title(['PICCS iter. ' num2str(outer) ' frame 115']); colormap gray; axis image; drawnow;
        subplot(2,2,2); plot(squeeze(abs(mean(errAll,3)))), title ('PICCS: Relative error norm respect to the full');
        s1=subplot(2,2,3:4); plot(squeeze(abs(u(22,20,:)))); title ('PICCS: Pixel Time course for pixel [22,20], the most significant');
        ylabel(s1,'Normalized units'),xlabel(s1,'Time frame');
        figure(h); waitbar(outer/nBreg,h);
        pause(1);
    end
    
end % outer

% Algorithm finished, update output
if (nargin >= 10)
    u = uBest;
end

% undo the normalization so that results are scaled properly
u           =   u/normFactor;
close(h);
return;


function normFactor = getNormalizationFactor(R,f)
normFactor = 1/norm(f(:)/size(R==1,1));
return;



function d = Dx(u)
[rows,cols,numTime] = size(u);
d = zeros(rows,cols,numTime);
d(:,2:cols,:) = u(:,2:cols,:)-u(:,1:cols-1,:);
d(:,1,:) = u(:,1,:)-u(:,cols,:);
return

function d = Dxt(u)
[rows,cols,numTime] = size(u);
d = zeros(rows,cols,numTime);
d(:,1:cols-1,:) = u(:,1:cols-1,:)-u(:,2:cols,:);
d(:,cols,:) = u(:,cols,:)-u(:,1,:);
return

function d = Dy(u)
[rows,cols,numTime] = size(u);
d = zeros(rows,cols,numTime);
d(2:rows,:,:) = u(2:rows,:,:)-u(1:rows-1,:,:);
d(1,:,:) = u(1,:,:)-u(rows,:,:);
return

function d = Dyt(u)
[rows,cols,numTime] = size(u);
d = zeros(rows,cols,numTime);
d(1:rows-1,:,:) = u(1:rows-1,:,:)-u(2:rows,:,:);
d(rows,:,:) = u(rows,:,:)-u(1,:,:);
return



function [xs,ys] = shrink2(x,y,lambda)

s = sqrt(x.*conj(x)+y.*conj(y));
ss = s-lambda;
ss = ss.*(ss>0);

s = s+(s<lambda);
ss = ss./s;

xs = ss.*x;
ys = ss.*y;

return;





