      program ex2
      implicit none
      integer e,maxX,m,l,A,B,i,n
      parameter(maxX=10,m=100,l=2,A=0,B=10)
      double precision delta,k,pi,Vp
      double precision, dimension(m+1) :: x,xi,V,w
      complex*16, dimension(m+1) :: Q,u
      complex j,kc,deltac,R,T,Rp,wc,c2
      parameter(j=(0,1))

      
      pi=dacos(-1d0)
      e=1
      n=m+1
      k=sqrt(dfloat(2*e))
      delta=dfloat(2*maxX)/dfloat(m)
      write(*,*)'delta = ',delta
      do i=1,n
        x(i)=-maxX+(i-1)*delta
        !write(*,*)'x (',i,') = ',x(i)
      enddo
      
      do i=1,n
        xi(i)=-exp(2d0*pi*x(i)/dfloat(l))
        !write(*,*)'xi (',i,') = ',xi(i)
        Vp=(dfloat(-A)*xi(i))/(1-xi(i))
        V(i)=Vp-dfloat(B)*xi(i)/((1-xi(i))**2d0)
        !write(*,*)'V (',i,') = ',V(i)
      enddo

      do i=1,n
        w(i)=2*(V(i)-dfloat(e))
        !write(*,*)'W (',i,') = ',w(i)
      enddo
      kc=CMPLX(k,0)
      deltac=CMPLX(delta,0)
      Q(n)=exp(j*kc*deltac)
      !write(*,*)'Q (n) = ',Q(n)
      c2=CMPLX(2,0)
      do i=n,2,-1
        wc=CMPLX(w(i),0)
        Q(i-1)=(wc*deltac**c2+c2-Q(i))**(CMPLX(-1,0))
        !write(*,*)'Q (',i-1,') = ',Q(i-1)
      enddo
      
      Rp=exp(j*kc*CMPLX(x(2),0))-Q(1)*exp(j*kc*CMPLX(x(1),0))

      R=Rp/(Q(1)*exp(-1*j*kc*CMPLX(x(1),0))-exp(-1*j*kc*CMPLX(x(2),0)))

      T=1-R*CONJG(R)
      write(*,*)'R = ',R
      u(1)=exp(-1*j*kc*CMPLX(x(1),0))+R*exp(j*kc*CMPLX(x(1),0))
      do i=1,n-1
        u(i+1)=Q(i)*u(i)
        write(*,*)'u (',i,') = ',u(i)
      enddo
      
      end
