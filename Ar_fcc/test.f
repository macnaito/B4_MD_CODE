      program test
      real*8 v(4116)
      open (11,file='final50.dat')
      do i=1,4116
       read(11,*)v(i)
      enddo
      close(11)
