
!! R.R -- spinor weight -- can check the phase of amplitudes
      spinorweight=1d0
      spinorweight(-1,-1,:,:)=spinorweight(-1,-1,:,:)*za(1,2)**2
      spinorweight(+1,+1,:,:)=spinorweight(+1,+1,:,:)*zb(1,2)**2
      spinorweight(-1,+1,:,:)=spinorweight(-1,+1,:,:)*za(1,3)**2*zb(3,2)**2
      spinorweight(+1,-1,:,:)=spinorweight(+1,-1,:,:)*zb(1,3)**2*za(3,2)**2

      spinorweight(:,:,-1,-1)=spinorweight(:,:,-1,-1)*za(3,5)*zb(4,6)
      spinorweight(:,:,-1,+1)=spinorweight(:,:,-1,+1)*za(3,6)*zb(4,5)
      spinorweight(:,:,+1,-1)=spinorweight(:,:,+1,-1)*za(4,5)*zb(3,6)
      spinorweight(:,:,+1,+1)=spinorweight(:,:,+1,+1)*za(4,6)*zb(3,5)
