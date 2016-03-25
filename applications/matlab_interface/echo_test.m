hippylib = HippyClient();

eta_len = 256;
eta=randn(eta_len, 1);

tic
for i = 1:2000
    hippylib.echo(eta)
end
toc
    
hippylib.close()