% Demo_Split_Bregman_PICCS_fMRI.m
%
% If you use this code, please, cite the following paper: 
% Chavarrias, C., Abascal, J.F., Montesinos, P. and Desco, M.,"Exploitation 
% of temporal redundancy in compressed sensing reconstruction of fMRI studies 
% with a prior-based algorithm (PICCS)". Med Phys,  42(7): p. 3814 (2015). 
% DOI:% http://dx.doi.org/10.1118/1.4921365
%
% Code downloaded from repository: 
% https://github.com/HGGM-LIM/Split_Bregman_PICCS_fMRI
% -------------------------------------------------------------------------
% Demo for reconstructing fMRI preclinical data with Prior Image-Based 
% Constrained Compressed Sensing (PICCS) using the Split Bregman formulation. 
% 
% PICCS is solved using a modification of Goldstein'n code mrics.m 
% downloaded from 
% (http://www.ece.rice.edu/~tag7/Tom_Goldstein/Split_Bregman.html), see Tom
% Goldstein and Stanley Osher. The Split Bregman Method for L1-Regularized
% Problems. SIAM J. Imaging Sci., 2(2), 323–343.  
%
% PICCS minimizes
% min_u (1-alpha)|grad_x,y u|_1 + alpha|grad_x,y(u-uprior)|_1 st. ||Fu-f||^2 < delta, proposed in 
% Chavarrias et al. Med Phys,  42(7): p. 3814 (2015).   
%
% -------------------------------------------------------------------------
% Data corresponds to the central slice of an fMRI study on a healthy rat 
% (originally 5 slices and 115 repetitions) with a block paradigm 
% OFF(15 volumes)-ON(5 volumes)-OFF-ON-OFF-ON-OFF-ON-OFF-ON-OFF. 
% Loaded data: Simulated absolute image from the detrended 4-coil kspaces 
% acquired 
% The acquired data set is available from
% http://biig.uc3m.es/XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%
% Undersampling is simulated using a modified version of Lustig's variable
% density pdf, downloaded from (SparseMRI V0.2)
% http://web.stanford.edu/~mlustig/SparseMRI.html, see M. Lustig, D.L
% Donoho and J.M Pauly "Sparse MRI: The Application of Compressed Sensing
% for Rapid MR Imaging" Magnetic Resonance in Medicine, 2007 Dec;
% 58(6):1182-1195.   
%
% Statistical analysis of the fMRI series reconstructed can be performed
% using the following software: https://github.com/HGGM-LIM/fmrat 
% If you use this software please cite the following article:
% Chavarrías, C., García-Vázquez, V., Alemán-Gómez, Y., Montesinos, P., 
% Pascau, J. and Desco, M.,
% "fMRat: an extension of SPM for a fully automatic analysis of rodent 
% brain functional magnetic resonance series". 
% Medical & Biological Engineering & Computing: p. 1-10 (2015)
%
% Cristina Chavarrías, Juan FPJ Abascal
% Departamento de Bioingeniería e Ingeniería Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% cristina.chavarrias@gmail.com, juanabascal78@gmail.com, desco@hggm.es



%==========================================================================
% LOAD DATA
%==========================================================================
cd(fileparts(which('Demo_Split_Bregman_PICCS_fMRI')));
load('fMRI_64x64x115_image.mat'); 

% Simulate data
data0       =   fft2(image0);         


%==========================================================================
% DISPLAY FULL DATASET, highest activation pixel is [22,20]
%==========================================================================
% Here we display the thresholded statistical map obtained from the 
% original fully sampled dataset analyzed with: https://github.com/HGGM-LIM/fmrat
% As well as the fully sampled 115 frames in the simulated fMRI series and 
% the time course of the pixel [22,20], which yielded the highest statistical 
% t value.

% Display the statistical map
    close('all');
    scrs    =   get(0,'ScreenSize');
    h1       =   figure; set(h1,'Position',[30,scrs(4)/2,scrs(3)/3,scrs(4)*0.8/2]);
    subplot(2,2,2); imshow('map_full_p_0_01_k_12.tif'); title('Objective: stats. map');

% Display the time frames    
    subplot(2,2,1); 
    for i=1:115 
        figure(h1);
        imagesc(abs(image0(:,:,i))); title(['Objective: image fr ' num2str(i) '/115']); 
        colormap gray; axis image; drawnow;
        pause(0.05);
    end
    
% Display the time course of the highest activation pixel    
    tcourse=subplot(2,2,3:4);  plot(squeeze(image0(22,20,:))); title(tcourse,'Objective: Pixel Time course for pixel [22,20], the most significant')
    ylabel(tcourse,'Absolute units'),xlabel(tcourse,'Time frame');

    
%==========================================================================
% UNDERSAMPLE DATASET, choose an undersampling rate (1/acc. factor)
%==========================================================================
%     sp                   = [ ...
%                         0.05;...      % 5% of lines   
%                         0.1; ...      % 10% of lines 
%                         0.125;...     % 12.5% of lines 
%                         0.15;...      % 15% of lines   
%                         0.175;...     % 17.5% of lines 
%                         0.2; ...      % 20% of lines 
%                         0.3; ...      % 30% of lines 
%                         0.4; ...      % 40% of lines 
%                         0.5; ...      % 50% of lines 
%                         0.6; ...      % 60% of lines 
%                         0.7; ...      % 70% of lines 
%                         0.8; ...      % 80% of lines 
%                         0.9 ...       % 90% of lines 
%                             ];
%    The following are parameters required to generate the pdf at these
%    undersampling rates above:
%     Ps          =   [17,10,8,8,7.5,5,3,2,2,2,2,2,2];        % Decay
%     N           =   size(image0,1);                         % #phase enc lines
%     radius      =   [2/N ,3/N, 3.3/N, 3.7/N, 4.1/N, ...     % Central radius
%                     4.5/N,7/N,7/N,7/N,7/N,7/N,7/N,7/N];
        
% To display the generated patterns for the 115 time frames
    h2          =   figure;
    set(h2,'Position',[30+scrs(3)/3,scrs(4)/2,scrs(3)/3,scrs(4)/2*0.8]);
    rname       =   'central_dense';
    
% Choose an undersampling rate (acc. factor)
    sp          =   0.3;     
    Ps          =   3;
    radius      =   7/size(image0,1);
    
% Create an undersampling for each time frame and display them
    dims        =   size(image0);
    RAll        =   genRAll(dims,sp,Ps,radius,pwd,rname,h2);
    fprintf('Generated RAll %s\t for sparsity: %s \r\n',rname, num2str(sp));                    

% Display the sum of the undersamplings along time to check the kspace completion 
    RAll_visu   =   fftshift(RAll);
    figure(h2),subplot(1,2,2),
    spy(squeeze(sum(RAll_visu,3))), 100*nnz(squeeze(sum(RAll_visu,3)))/size(image0,1)
    title('Sum of all patterns accross time -to check completion-');

    
%==========================================================================
% Apply undersampling
%==========================================================================
    data        = data0.*RAll;
    % ---------------------------------------------------------------------
    % IFFT, the zero filling reference
    u_ifft      = ifft2(data);

    
%==========================================================================
% Prepair algorithm inputs and output
%==========================================================================
    
%reco pars
    dimIm               =   size(image0);
    muu                 =   1;      % fidelity term
    lambda              =   1;
    alpha               =   0.95;    % Prior weight // (1-alpha) for Dxy(u)
    gamma               =   1;      % +gamma*u  // stability term
    nBreg               =   150;
    im_cmp              =   zeros(size(image0));
    errAll              =   zeros(1,nBreg);


% Build prior as the temporal mean of all undersampled kspaces
    av                          =   squeeze(sum(data,3));
    idx                         =   sum(RAll,3);
    prior                       =   zeros(size(av));
    prior(find(idx(:,1)~=0),:)  =   av(find(idx(:,1)~=0),:)./idx(find(idx(:,1)~=0),:);  % normalize the lines intensities according to averaged lines
    uprior                      =   ifft2(prior);
    uprior                      =   repmat(uprior,[1,dimIm(3)]);
    uprior                      =   reshape(uprior,dimIm); 
    
    RAll                        =   (data~=0);   % avoid zeros in the full data, e.g. EPI boundaries

 % For visualization of the reconstructed image evolution with iterations
    h4                          =   figure; 
    set(h4,'Position',[30,30,scrs(3)/3,scrs(4)/2*0.8]);    
    
    
%==========================================================================
% PICCS CALL
%==========================================================================
    
% Reconstruction and display PICCS results at each iteration
    tic
    [u_PICCS, err_PICCS]    =   PICCS_fMRI(RAll, data, dimIm, uprior, muu, lambda,alpha, gamma, nBreg, image0, h4);
    toc

    figure(h4), 
    tcourse                 =   subplot(2,2,3:4); 
    plot(squeeze(abs(u_PICCS(22,20,:)))); 
    title (tcourse,'PICCS: Pixel Time course for pixel [22,20], the most significant');
    ylabel(tcourse,'Absolute units');xlabel(tcourse,'Time frame')

% Show image difference (ON-OFF). For an accurate activation map you
    % must perform a full analysis (corrections, statistics, etc) with SPM, 
    % FSL or the preclinical tool fMRat: https://github.com/HGGM-LIM/fmrat
    stim                =   zeros([115,1]);
    stim([16:20,36:40,56:60,76:80,96:100],1)    =   1;
    diff                =   abs(mean(u_PICCS(:,:,squeeze(logical(stim))),3)-mean(u_PICCS(:,:,squeeze(logical(~stim))),3));
    h5                      =   figure;
    set(h5,'Position',[30+scrs(3)/6,30+scrs(4)/4*0.8,scrs(3)/6,scrs(4)/4*0.8]);
    imagesc(diff); colormap gray; axis image; drawnow;
    title ('PICCS: Difference (ON-OFF)');    
    
    
%==========================================================================    
% Show ZERO FILLING results
%==========================================================================    
    
% Show frame 115
    h6                      =   figure;
    set(h6,'Position',[30+scrs(3)/3,30,scrs(3)/3,scrs(4)/2*0.8]);
    subplot(2,2,1); 
    imagesc(squeeze(abs(u_ifft(:,:,115)))); colormap gray; axis image; drawnow;
    title ('Zero filled: Frame 115');
    
% Show the [22,20] pixel time course
    tcourse                 =    subplot(2,2,3:4); 
    plot(squeeze(abs(u_ifft(22,20,:)))); 
    title (tcourse,'Zero filled: Pixel Time course for pixel [22,20], the most significant');
    ylabel(tcourse,'Absolute units');xlabel(tcourse,'Time frame')
    
% Show the difference ON-OFF
    diff                =   abs(mean(u_ifft(:,:,squeeze(logical(stim))),3)-mean(u_ifft(:,:,squeeze(logical(~stim))),3));
    subplot(2,2,2); 
    imagesc(diff); colormap gray; axis image; drawnow;
    title ('Zero filled: Difference (ON-OFF)');        
