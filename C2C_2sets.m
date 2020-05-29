function matPred = C2C_2sets(conMat2useG1, conMat2predG1, conMat2useG2, nCompPCA, nCompPLS, meanPCA)
% conMat : subjects x edges

if nargin < 4;    nCompPCA = [];    end;    
if nargin < 5;    nCompPLS = [];    end;    
if nargin < 6;    meanPCA = 0;    end;

% nCompPCA=100;    nCompPLS = 10;    meanPCA = 0;

%% Run C2C pipeline
disp('    Running C2C transformation');

% ----  Construct a state transition model using a training set of group 1

    % Define subnetwork templates for each state in training set with PCA
    if meanPCA == 0
        [s1coeff, s1score, ~, ~, ~, s1mu] = pca(conMat2useG1, 'Centered', false);        [s2coeff, s2score, ~, ~, ~, s2mu] = pca(conMat2predG1, 'Centered', false);
    elseif meanPCA == 1
        [s1coeff, s1score, ~, ~, ~, s1mu] = pca(conMat2useG1);        [s2coeff, s2score, ~, ~, ~, s2mu] = pca(conMat2predG1);
    end;    

    % Select the number of top subnetwork to use for state transition
    if ~isempty(nCompPCA)
        s1score(:, [(nCompPCA+1):end]) = [];    s2score(:, [(nCompPCA+1):end]) = [];
        s1coeff(:, [(nCompPCA+1):end]) = [];    s2coeff(:, [(nCompPCA+1):end]) = [];
    end
    
    % Define a function of state transformation with PLSR of specific number of PLSR components
    if ~isempty(nCompPLS);        [~,~,~,~, betaPLS] = plsregress( s1score, s2score, nCompPLS);
    else;                                   [~,~,~,~, betaPLS] = plsregress( s1score, s2score);
    end;

% ----  Apply the constructed transition model to a testing set of group 2

    % Estimate, in testing set, weights for the predefined subnetwork templates of state 1
    if meanPCA == 0;              s1score2pred = conMat2useG2 / s1coeff';
    elseif meanPCA == 1;        s1score2pred = (conMat2useG2 - s1mu) / s1coeff';    end;

    % Predict, in testing set, weights for the predefined subnetwork templates of state 2
    s2scorePred = [ones(size(s1score2pred,1),1) s1score2pred]  * betaPLS;    clear beta

    % Predict, in testing set, the whole brain connectome of state 2
    if meanPCA == 0;              matPred = (s2scorePred * s2coeff');
    elseif meanPCA == 1;        matPred = (s2scorePred * s2coeff') + s2mu;    
    end;    clear s2scorePred
    
disp('    ');

