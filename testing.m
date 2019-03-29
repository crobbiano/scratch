load('frm_trex_pond_insitu_2class_test_insitu_single_gen_ol_0.1.mat')

%use test0data
D=test0data(:,1:1000);
D=Unitary(271);
G=test0data(:,1001:end);
V=D\G;
covG = cov(G');

[T,~]= eig(covG);
VV=OMP(T'*D,G,15,.01);