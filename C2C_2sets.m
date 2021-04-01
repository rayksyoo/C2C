function matPred = C2C_2sets(conMat2useG1, conMat2predG1, conMat2useG2, nCompPCA, nCompPLS, meanPCA)
% conMat : subjects x edges (the number of subejcts by the number of 1-dim vectorized connectivity edges)

if nargin < 4;    nCompPCA = 100;    end;    
if nargin < 5;    nCompPLS = 10;    end;    
if nargin < 6;    meanPCA = 0;    end;

%% Run C2C pipeline
disp('    Running C2C transformation ...');

% ----  Construct a state transformation model in a training set (Group 1)

    % Define subnetworks for each state (resting state and task-related state) in training set with PCA
    if meanPCA == 0
        [s1coeff, s1score, ~, ~, ~, s1mu] = pca(conMat2useG1, 'Centered', false);        [s2coeff, s2score, ~, ~, ~, s2mu] = pca(conMat2predG1, 'Centered', false);
    elseif meanPCA == 1
        [s1coeff, s1score, ~, ~, ~, s1mu] = pca(conMat2useG1);        [s2coeff, s2score, ~, ~, ~, s2mu] = pca(conMat2predG1);
    end;    

    % Select the number (nCompPCA) of top subnetworks to use for state transformation
    if ~isempty(nCompPCA)
        s1score(:, [(nCompPCA+1):end]) = [];    s2score(:, [(nCompPCA+1):end]) = [];
        s1coeff(:, [(nCompPCA+1):end]) = [];    s2coeff(:, [(nCompPCA+1):end]) = [];
    end
    
    % Define a state transformation with PLSR with a specific number of PLSR components (nCompPLS)
    if ~isempty(nCompPLS);        [~,~,~,~, betaPLS] = plsregress( s1score, s2score, nCompPLS);
    else;                                   [~,~,~,~, betaPLS] = plsregress( s1score, s2score);
    end;

% ----  Apply the constructed state transformation model to a testing set (Group 2)

    % Estimate weights for the predefined subnetworks of state 1 (e.g., resting-state subnetworks) of a testing set
    if meanPCA == 0;              s1score2pred = conMat2useG2 / s1coeff';
    elseif meanPCA == 1;        s1score2pred = (conMat2useG2 - s1mu) / s1coeff';    end;

    % Predict weights for the predefined subnetworks of state 2 (e.g., task-related subnetworks) of a testing set
    s2scorePred = [ones(size(s1score2pred,1),1) s1score2pred]  * betaPLS;

    % Predict the whole brain connectome of state 2 of a testing set
    if meanPCA == 0;              matPred = (s2scorePred * s2coeff');
    elseif meanPCA == 1;        matPred = (s2scorePred * s2coeff') + s2mu;    
    end
    
disp('    Connectomes generated.');
