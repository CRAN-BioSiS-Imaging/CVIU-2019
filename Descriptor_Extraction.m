function vectorF = Descriptor_Extraction(I, para)
Descriptor = para.D;
[m,n] = size(I);

switch Descriptor
%% BINARY DESCRIPTORS (BD):
    %% Census transform:
    case 'BD_Census'
        % R. Zabih, J. Woodfill, 
        % Non-parametric local transforms for computing visual correspondence, ECCV 1994, pp. 151–158
        vectorF = zeros(m,n,8);
        B = zeros(m-2,n-2,9);
        k=0;
        for i = 1:3
            for j = 1:3
                k = k + 1;
                B(:,:,k) = I(i:(i+m-3),j:(j+n-3));
            end
        end
        k = 0;
        for t = 1:9
            if t ~= 5
                k = k+1;
                vectorF(2:end-1,2:end-1,k) = B(:,:,t) - B(:,:,5) ;
            end
        end
        vectorF = (vectorF > 0);
    
    %% CRT transform:
    % O. Demetz, D. Hafner, J. Weickert, 
    % "The complete rank transform: A tool for accurate and morphologically
    % invariant matching of structures,"
    % BMVC 2013.
    case 'BD_CRT' 
        vectorF = zeros(m,n,9);
        B = zeros(m-2,n-2,9);
        k=0;
        for i = 1:3
            for j = 1:3
                k = k + 1;
                B(:,:,k) = I(i:(i+m-3),j:(j+n-3));
            end
        end
        
        for t = 1:9
            M = zeros(m-2,n-2);
            for tt = 1:9
                if tt ~=t
                    D6 = (B(:,:,t) > B(:,:,tt)) ;
                    M = M + D6;
                end
            end
            vectorF(2:end-1,2:end-1,t) = M;
        end
        
    %% LDP transform:
    %M. H. Kabir, T. Jabid, O. Chae, 
    % A local directional pattern variance (LDPv) based face descriptor for human facial expression recognition,
    % IEEE International Conference on Advanced Video and Signal Based Surveillance, 2010, pp. 526–532
    case 'D_LDP'
        vectorF = zeros(m,n,8);
        M1 = [-3 -3 5;-3 0 5; -3 -3 5]; M2 = [-3 5 5; -3 0 5; -3 -3 -3];
        M3 = [5 5 5; -3 0 -3; -3 -3 -3]; M4 = [5 5 -3; 5 0 -3; -3 -3 -3];
        M5 = [5 -3 -3; 5 0 -3; 5 -3 -3]; M6 = [-3 -3 -3; 5 0 -3; 5 5 -3];
        M7 = [-3 -3 -3;-3 0 -3; 5 5 5]; M8 = [-3 -3 -3; -3 0 5; -3 5 5];
        F(:,:,1) = imfilter(I,M1);
        F(:,:,2) = imfilter(I,M2);
        F(:,:,3) = imfilter(I,M3);
        F(:,:,4) = imfilter(I,M4);
        F(:,:,5) = imfilter(I,M5);
        F(:,:,6) = imfilter(I,M6);
        F(:,:,7) = imfilter(I,M7);
        F(:,:,8) = imfilter(I,M8);
        
        Fabs = abs(F);
        Fsort = sort(Fabs,3);
        k_largest = Fsort(:,:,3); % k=3
        
        for k = 1:8
            vectorF(:,:,k) = (Fabs(:,:,k) > k_largest);
        end
        
    %% MLDP transform:
    % M. A. Mohamed, H. A. Rashwan, B. Mertsching, M. A. Garc´?a, D. Puig,
    % "Illumination-robust optical flow using a local directional pattern," 
    % IEEE Transactions on Circuits and Systems for Video Technology 24 (9) (2014)1499–1508
    case 'BD_MLDP'
        M1 = [-3 -3 5;-3 0 5; -3 -3 5]; M2 = [-3 5 5; -3 0 5; -3 -3 -3];
        M3 = [5 5 5; -3 0 -3; -3 -3 -3]; M4 = [5 5 -3; 5 0 -3; -3 -3 -3];
        M5 = [5 -3 -3; 5 0 -3; 5 -3 -3]; M6 = [-3 -3 -3; 5 0 -3; 5 5 -3];
        M7 = [-3 -3 -3;-3 0 -3; 5 5 5]; M8 = [-3 -3 -3; -3 0 5; -3 5 5];
        F(:,:,1) = imfilter(I,M1);
        F(:,:,2) = imfilter(I,M2);
        F(:,:,3) = imfilter(I,M3);
        F(:,:,4) = imfilter(I,M4);
        F(:,:,5) = imfilter(I,M5);
        F(:,:,6) = imfilter(I,M6);
        F(:,:,7) = imfilter(I,M7);
        F(:,:,8) = imfilter(I,M8);
        vectorF = sign(F);
    
    %% Proposed descriptor:
    case 'BD_new1' %% Robinson mask
        % Dinh-Hoan Trinh, Christian Daul,
        % "On Illumination-Invariant Variational Optical Flow for Weakly
        % Textured Scenes," CVIU 2019.
        vectorF = zeros(m,n,8);
        M1 = [-1 0 1;-2 0 2; -1 0 1]; M2 = [0 1 2; -1 0 1; -2 -1 0];
        M3 = [1 2 1; 0 0 0; -1 -2 -1]; M4 = [2 1 0; 1 0 -1; 0 -1 -2];
        M5 = [1 0 -1; 2 0 -2; 1 0 -1]; M6 = [0 -1 -2; 1 0 -1; 2 1 0];
        M7 = [-1 -2 -1;0 0 0; 1 2 1]; M8 = [-2 -1 0; -1 0 1; 0 1 2];
        F(:,:,1) = imfilter(I,M1);
        F(:,:,2) = imfilter(I,M2);
        F(:,:,3) = imfilter(I,M3);
        F(:,:,4) = imfilter(I,M4);
        F(:,:,5) = imfilter(I,M5);
        F(:,:,6) = imfilter(I,M6);
        F(:,:,7) = imfilter(I,M7);
        F(:,:,8) = imfilter(I,M8);
        vectorF = sign(F);

%% REAL VALUE DESCRIPTORS (RD):
    %% 1. Correlation transform:
    case 'RD_Corr'
        %M. Drulea, S. Nedevschi, 
        %"Motion estimation using the correlation transform," 
        %IEEE Transactions on Image Processing 22 (8) (2013) 3260–3270.
        vectorF = zeros(m,n,9);
        B = zeros(m-2,n-2,9);
        k=0;
        for i = 1:3
            for j = 1:3
                k = k + 1;
                B(:,:,k) = I(i:(i+m-3),j:(j+n-3));
            end
        end
        
        Mean = mean(B,3);
        Std  = std(B,1,3) + 1e-09;
        for k = 1:9
            vectorF(2:end-1,2:end-1,k) = (B(:,:,k) - Mean)./Std;
        end
        
    %% 2.NND transform:
    case 'RD_NND'
        %S. Ali, C. Daul, E. Galbrun, W. Blondel:
        %"Illumination invariant optical flow using neighborhood descriptors". 
        %Computer Vision and Image Understanding 145: 95-110 (2016)
        vectorF = RF_NND_feature_vector(I,para);
    
    %% 3. New1
    case 'RD_IIOF'%% Robinson mask
        %Dinh Hoan Trinh, Walter Blondel, Christian Daul:
        %A general form of illumination-invariant descriptors in variational optical flow estimation. 
        %ICIP 2017: 2533-2537
        vectorF = zeros(m,n,8);
        M1 = [-1 0 1;-2 0 2; -1 0 1]; M2 = [0 1 2; -1 0 1; -2 -1 0];
        M3 = [1 2 1; 0 0 0; -1 -2 -1]; M4 = [2 1 0; 1 0 -1; 0 -1 -2];
        M5 = [1 0 -1; 2 0 -2; 1 0 -1]; M6 = [0 -1 -2; 1 0 -1; 2 1 0];
        M7 = [-1 -2 -1;0 0 0; 1 2 1]; M8 = [-2 -1 0; -1 0 1; 0 1 2];
        F(:,:,1) = imfilter(I,M1);
        F(:,:,2) = imfilter(I,M2);
        F(:,:,3) = imfilter(I,M3);
        F(:,:,4) = imfilter(I,M4);
        F(:,:,5) = imfilter(I,M5);
        F(:,:,6) = imfilter(I,M6);
        F(:,:,7) = imfilter(I,M7);
        F(:,:,8) = imfilter(I,M8);
        Nor = F(:,:,1).^2 + F(:,:,2).^2 + F(:,:,3).^2 + F(:,:,4).^2 + F(:,:,5).^2 + F(:,:,6).^2 + F(:,:,7).^2 + F(:,:,8).^2;
        Nor = sqrt(Nor) + 1e-09;
        for k = 1:8
            vectorF(:,:,k) = F(:,:,k)./Nor;
        end
        
    %% 4. New2
    % Dinh-Hoan Trinh, Christian Daul,
    % "On Illumination-Invariant Variational Optical Flow for Weakly
    % Textured Scenes," CVIU 2019.
    case 'RD_new2'
        vectorF = zeros(m,n,9);
        k=0;
        B = zeros(m-2,n-2,9);
        for i = 1:3
            for j = 1:3
                k = k + 1;
                B(:,:,k) = I(i:(i+m-3),j:(j+n-3));
            end
        end
        Vmax = max(B,[],3);
        Vmin = min(B,[],3);
        
        Vmax_Vmin = (Vmax - Vmin) +  1e-09;
        for k = 1:9
            B(:,:,k) = (B(:,:,k) - Vmin)./Vmax_Vmin ;
        end
        vectorF(2:end-1,2:end-1,:)  = exp(B);
    
    %% 5. New 3:
    % D.-H. Trinh, W. Blondel, D. Lamarque, C. Daul, 
    % "Illumination-invariant optical flow: Application to endoscopic image mosaicing," 
    % in: XXVIe Colloque GRETSI Traitement du Signal et des Images, GRETSI 2017.
    case 'RD_Gretsi'
        vectorF = zeros(m,n,8);
        M1 = [-3 -3 5;-3 0 5; -3 -3 5]; M2 = [-3 5 5; -3 0 5; -3 -3 -3];
        M3 = [5 5 5; -3 0 -3; -3 -3 -3]; M4 = [5 5 -3; 5 0 -3; -3 -3 -3];
        M5 = [5 -3 -3; 5 0 -3; 5 -3 -3]; M6 = [-3 -3 -3; 5 0 -3; 5 5 -3];
        M7 = [-3 -3 -3;-3 0 -3; 5 5 5]; M8 = [-3 -3 -3; -3 0 5; -3 5 5];
        F(:,:,1) = imfilter(I,M1);
        F(:,:,2) = imfilter(I,M2);
        F(:,:,3) = imfilter(I,M3);
        F(:,:,4) = imfilter(I,M4);
        F(:,:,5) = imfilter(I,M5);
        F(:,:,6) = imfilter(I,M6);
        F(:,:,7) = imfilter(I,M7);
        F(:,:,8) = imfilter(I,M8);
        Nor = F(:,:,1).^2 + F(:,:,2).^2 + F(:,:,3).^2 + F(:,:,4).^2 + F(:,:,5).^2 + F(:,:,6).^2 + F(:,:,7).^2 + F(:,:,8).^2;
        Nor = sqrt(Nor) + 1e-09;
        for k = 1:8
            vectorF(:,:,k) = F(:,:,k)./Nor;
        end
        
    %% 6. New 4    
    case 'RD_PR2020'
        vectorF = zeros(m,n,12);
        M1 = [-1 -1 -1; 0 3 0; 0 0 0]; F(:,:,1) = imfilter(I,M1);
        M2 = [0 -1 -1; 0 3 -1; 0 0 0]; F(:,:,2) = imfilter(I,M2);
        M3 = [0 0 -1; 0 3 -1; 0 0 -1]; F(:,:,3) = imfilter(I,M3);
        M4 = [0 0 -1; 0 3 -1; 0 -1 -1]; F(:,:,4) = imfilter(I,M4);
        M5 = [0 0 0; 0 3 0; -1 -1 -1]; F(:,:,5) = imfilter(I,M5);
        M6 = [0 0 0; -1 3 0; -1 -1 0]; F(:,:,6) = imfilter(I,M6);
        M7 = [-1 0 0; -1 3 0; -1 0 0]; F(:,:,7) = imfilter(I,M7);
        M8 = [-1 -1 0; -1 3 0; 0 0 0]; F(:,:,8) = imfilter(I,M8);
        M9 = [0 -1 0; -1 3 -1; 0 0 0]; F(:,:,9) = imfilter(I,M9);
        M10 = [0 -1 0; 0 3 -1; 0 -1 0]; F(:,:,11) = imfilter(I,M10);
        M11 = [0 0 0; -1 3 -1; 0 -1 0]; F(:,:,12) = imfilter(I,M11);
        M12 = [0 -1 0; -1 3 0; 0 -1 0]; F(:,:,13) = imfilter(I,M12);
        
        Norm_F = zeros(m,n);
        for i = 1:12
            Norm_F = Norm_F + F(:,:,i).^2;
        end
        Norm_F = sqrt(Norm_F) + 1e-09;
        
        for k = 1:12
            vectorF(:,:,k) = F(:,:,k)./Norm_F;
        end
        
        
        
    otherwise
        disp('Do not know Descriptor!');
end






