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
disp(negLogPost);
hippylib.close();