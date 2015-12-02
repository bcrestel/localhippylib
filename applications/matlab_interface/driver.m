close all
clear all
clear classes

hippylib = HippyClient();
hippylib.computeMapPoint();
disp('Computed map point');
kdim = hippylib.KLE_GaussianPost();
disp('Computed KLE Gaussian Posterior');
hippylib.close();