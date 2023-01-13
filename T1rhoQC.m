%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Authors: Deep B. Gandhi and Jonathan Dudley
%Last modified date: 5/23/2019
%Last modified by: Deep B. Gandhi 
%Purpose: To perform QC over T1rho images
%Inputs: T1rho images with different spin lock time
%Outputs: A png image file with QC report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subID = ''; % String indicating subject ID


errmsg = '';
trace = cell(4,1);
for sliceN = 1:4
    % A 1x6 cell array of strings providing the full file names of the six
    % spin lock time images for a given slice
    TSLs = {'','','','','',''};
    
    % Check that slice prescriptions do match
    heads = cellfun(@dicominfo,TSLs);
    delPos = diff([heads.ImagePositionPatient],1,2);
    delOr = diff([heads.ImageOrientationPatient],1,2);
    if any(delPos(:)) || any(delOr(:))
        fprintf(2,'Slice prescriptions do not match across spin lock times!\n');
        errmsg = ': Slice prescriptions do not match across spin lock times!';
    end

    trace_threshold1 = .25;
    trace_threshold2 = .5;

    A = mat2gray(dicomread(TSLs{1}));
    B1 = bwperim(mat2gray(dicomread(TSLs{2}))>trace_threshold1);
    C1 = bwperim(mat2gray(dicomread(TSLs{3}))>trace_threshold1);
    D1 = bwperim(mat2gray(dicomread(TSLs{4}))>trace_threshold1);
    E1 = bwperim(mat2gray(dicomread(TSLs{5}))>trace_threshold1);
    F1 = bwperim(mat2gray(dicomread(TSLs{6}))>trace_threshold1);
    B2 = bwperim(mat2gray(dicomread(TSLs{2}))>trace_threshold2);
    C2 = bwperim(mat2gray(dicomread(TSLs{3}))>trace_threshold2);
    D2 = bwperim(mat2gray(dicomread(TSLs{4}))>trace_threshold2);
    E2 = bwperim(mat2gray(dicomread(TSLs{5}))>trace_threshold2);
    F2 = bwperim(mat2gray(dicomread(TSLs{6}))>trace_threshold2);
    trace{sliceN} = gray2ind(cat(2,...
        cat(3,A+B1/2+B2,A+B1/2,A),...
        cat(3,A+C1/2+C2,A+C1/2,A),...
        cat(3,A+D1/2+D2,A+D1/2,A),...
        cat(3,A+E1/2+E2,A+E1/2,A),...
        cat(3,A+F1/2+F2,A+F1/2,A)),...
        256);
end

imshow(cat(1,trace{1},trace{2},trace{3},trace{4}))
title([subID errmsg])
xlabel('Spin Lock Times -->','VerticalAlignment','Bottom')
ylabel('<-- Slices -->','VerticalAlignment','Top')
imwrite(frame2im(getframe(gcf)),[subID '_QC.png']);