
[ze,pe,me] = spinopidx(clusterSize,iSpin);
[ez,ep,em] = spinopidx(clusterSize,jSpin);
[zz,rl,lr,zr,zl,rz,lz,rr,ll] = spinopidx2(clusterSize,iSpin,jSpin);
doTest = true;
max(max(abs( ...
  SpinOp(:,:,z) - 2*commutator(SpinOp(:,:,rz),SpinOp(:,:,lz) ) )))
