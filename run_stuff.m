% Load array of OMP-MSC thresholded statistics, i.e., class labels, order
% sequentially
labels = load(labels.mat);

%% fake stuff
labels = [1, 2, 2, 1, 4, 1, 1];
% Initialize the Dirichlet-Categorical model
numClasses = 5;
alphas = ones(numClasses, 1);
%%
% Loop over the labels
ig = zeros(length(labels),1);
for i = 1:length(labels)
    % calculate the IG 
    [ig(i), alphas] = calcIG(alphas, labels(i));
end
