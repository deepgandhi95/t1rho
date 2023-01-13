%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Authors: Deep B. Gandhi
%Last modified date: 10/24/2019
%Last modified by: Deep B. Gandhi
%Purpose: To perform T1 rho parametric mapping
%Inputs: T1rho images with different spin-lock times
%Outputs: T1rho parametric map 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directory = pwd;
subid = dir(fullfile(directory));

for h = 1:numel(subid)
    T1rho_filepath = convertCharsToStrings(fullfile(subid(h).folder,subid(h).name));
    cd(T1rho_filepath); 
    subID = subid(h).name;
    File_path_New = strcat(T1rho_filepath,'\NewPro_Niftifiles');
    
    if exist(File_path_New,'dir')
    
        I = dir('NewPro_Niftifiles');
        I(ismember( {I.name}, {'.', '..'})) = [];
        Niftifiles = convertCharsToStrings(I(1).folder);
        cd(Niftifiles)
        data = dir('*TSL*.nii');
        N = natsortfiles({data.name});
        TestImage1 = niftiread('Liver_regTSL0_Mask.nii'); 
        [x,y] = size(TestImage1(:,:,1)); % number of pixels in x and y

        % Read all the nifti files and concatenate into a single matrix
        t1rho_image = zeros(x,y,length(N),size(TestImage1,3)); % zeros in 4D 
        for n=1:length(N) %read all the nifti files
           t1rho_image(:,:,n,:) = niftiread(N{n});
           if size(t1rho_image(:,:,n,2)) ~= size(t1rho_image(:,:,n,1))
               disp('error: image matrix input sizes not equal');
           else 
                % do nothing 
           end 
        end


        % input the Spin-lock times
        tsl = [0 10 20 30 40 80];


        for p = 1:size(t1rho_image(:,:,:,:),4)
            Img = t1rho_image(:,:,:,p);
            [x,y,nEcho] = size(Img); % dimension of the Image  
            ymat = reshape(Img,[],nEcho);
            Ydata = reshape(ymat',[],1);
            % xdata is a vector which contains the echo time for each pixel element
            % in Y
            xmat = kron(tsl,ones(x*y,1));
            Xdata = reshape(xmat',[],1); 
            % It is necessary to index our problem. i.e. in our 1D array which
            % element corresponds to which problem.
            batchindex = kron(1:x*y,ones(1,nEcho)); 

            estcoefs1 = zeros(x*y,1);
            estcoefs2 = zeros(x*y,1);% estimated coefficients a & b
            % we don't want to go all the way through matrix. It is sufficient to
            % say that if a pixel is blank at each echo time then it will not have a
            % T1rho value.
            msk = sum(Img,3);
            mask = reshape(msk,[],1)'>0;

            parfor i = 1:x*y
                if  mask(i) == 1 % checks the echo times have non-zero pixel values
                    k = (batchindex == i);
                    Xk = Xdata(k); % this is the Echo time data for problem i 
                    Yk = Ydata(k); % this is the pixel value data for problem i
                    Yk = Yk(:)/max(Yk(:)); % normalise the Y-data for each series
                    Deep = find(Yk == 0);
                    if Deep~=0
                        Yk(Deep) = [];
                        Xk(Deep) = [];
                    else 
                    end

                    [f2,gof,output] = fit(Xk,Yk,'exp1');
                    T1_rho = -1/f2.b;
                    AdjRsquare = gof.adjrsquare;
                    estcoefs1(i) = T1_rho;
                    estcoefs2(i) = AdjRsquare;
                end
            end
            T1_rho1 = estcoefs1;
            AdjRsquare = estcoefs2;
            T1rhoimg(:,:,p) = reshape(T1_rho1',x,y); 
            AdjRsquareimg(:,:,p) = reshape(AdjRsquare',x,y);

        end
        niftiwrite(T1rhoimg,'T1rhomap.nii');
        niftiwrite(AdjRsquareimg,'AdjRsquaremap.nii');
        Namastedeep = zeros(size(T1rhoimg));
        for i = 1:size(T1rhoimg,3)
            figure; I1 = imshow(T1rhoimg(:,:,i),[]);
            rhoval = T1rhoimg < 0.45;
            rhoads = AdjRsquareimg < 0.98;
            Namaste = ~(rhoval | rhoads);
            T1rhovalue(i) = mean(nonzeros(T1rhoimg(:,:,i).*Namaste(:,:,i)));
            T1rhostddev(i) = std(nonzeros(T1rhoimg(:,:,i).*Namaste(:,:,i)));
        end
        save('T1rhomean.mat','T1rhovalue');
        save('T1rhostddev.mat','T1rhostddev');
        rehash toolboxcache
    end
end