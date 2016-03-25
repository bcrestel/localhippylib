close all
clear all
clear classes

hippylib = HippyClient();
hippylib.computeMapPoint();
disp('Computed map point');
kdim = hippylib.KLE_GaussianPost();
disp('Computed KLE Gaussian Posterior');
eta = zeros(kdim,1);
negLogPost = hippylib.negLogPost(eta);
negLogGaussianPost = hippylib.negLogGaussianPost(eta);
negLogLikelihood = hippylib.negLogLikelihood(eta);
negLogPrior = hippylib.negLogPrior(eta);
fprintf('negLogPost = %f, negLogLike = %f, negLogPrior = %f\n', negLogPost, negLogLikelihood, negLogPrior);
fprintf('negLogGaussianPost = %f\n', negLogGaussianPost);
hippylib.close();