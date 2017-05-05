n=100

D0=rep(-1,n);
D1=rep(1,n-1);
D0[1]=D0[n]=1;
D1[1]=D1[n-1]=0;
A=cbind(c(D1,0),D0,c(D1,0));

D=-1*diag(n)
D[2:n,1:(n-1)]=diag(n-1);
D[1,1]=D[n,n]=1
D[1,2]=D[2,1]=D[n,n-1]=D[n-1,n]=0
test_that("Rotating matices...", {
  expect_equal(rot2mat(A),D)
})