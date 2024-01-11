  MODULE parcel_module

  implicit none

  private
  public :: parcel_driver,parcel_interp,parcel_write,setup_parcel_vars

  CONTAINS

      subroutine parcel_driver(dt,xh,uh,ruh,xf,yh,vh,rvh,yf,zh,mh,rmh,zf,mf,    &
                               znt,rho,ua,va,wa,pdata)
      use input
      use constants
      use bc_module
      use comm_module
#ifdef MPI
      use mpi
#endif
      implicit none

!-----------------------------------------------------------------------
!  This subroutine updates the parcel locations
!-----------------------------------------------------------------------

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je) :: znt
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(nparcels,npvals) :: pdata

      integer :: n,np,i,j,k,iflag,jflag,kflag
      real :: uval,vval,wval,rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8,wsum
      real :: rxu,ryv,rzw,rxs,rys,rzs
      real :: x3d,y3d,z3d
      integer :: nrkp
      real :: dt2,uu1,vv1,ww1
      real :: z0,rznt,var

      logical, parameter :: debug = .false.

!----------------------------------------------------------------------
!  apply bottom/top boundary conditions:
!  [Note:  for u,v the array index (i,j,0) means the surface, ie z=0]
!     (for the parcel subroutines only!)

!$omp parallel do default(shared)  &
!$omp private(i,j)
  DO j=jb,je+1

    IF(bbc.eq.1)THEN
      ! free slip ... extrapolate:
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,0) = cgs1*ua(i,j,1)+cgs2*ua(i,j,2)+cgs3*ua(i,j,3)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,0) = cgs1*va(i,j,1)+cgs2*va(i,j,2)+cgs3*va(i,j,3)
      enddo
    ELSEIF(bbc.eq.2)THEN
      ! no slip:
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,0) = 0.0
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,0) = 0.0
      enddo
    ELSEIF(bbc.eq.3)THEN
      ! u,v near sfc are determined below using log-layer equations
    ENDIF

!----------

    IF(tbc.eq.1)THEN
      ! free slip ... extrapolate:
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,nk+1) = cgt1*ua(i,j,nk)+cgt2*ua(i,j,nk-1)+cgt3*ua(i,j,nk-2)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,nk+1) = cgt1*va(i,j,nk)+cgt2*va(i,j,nk-1)+cgt3*va(i,j,nk-2)
      enddo
    ELSEIF(tbc.eq.2)THEN
      ! no slip:
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,nk+1) = 0.0
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,nk+1) = 0.0
      enddo
    ENDIF

!----------

      ! assuming no terrain:
      IF(j.le.je)THEN
      do i=ib,ie
        wa(i,j,1)    = 0.0
        wa(i,j,nk+1) = 0.0
      enddo
      ENDIF

  ENDDO

!----------------------------------------------------------------------
!  Loop through all parcels:  if you have it, update it's location:

    dt2 = dt/2.0

    nploop:  &
    DO np=1,nparcels

      x3d = pdata(np,prx)
      y3d = pdata(np,pry)
      z3d = pdata(np,prz)

      iflag = -100
      jflag = -100
      kflag = 0

  ! cm1r19:  skip if we already know this processor doesnt have this parcel
  haveit1:  &
  IF( x3d.ge.xf(1) .and. x3d.le.xf(ni+1) .and.  &
      y3d.ge.yf(1) .and. y3d.le.yf(nj+1) )THEN

    IF(nx.eq.1)THEN
      iflag = 1
    ELSE
      ! cm1r19:
      i = ni+1
      do while( iflag.lt.0 .and. i.gt.1 )
        i = i-1
        if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
          iflag = i
        endif
      enddo
    ENDIF

    IF(axisymm.eq.1.or.ny.eq.1)THEN
      jflag = 1
    ELSE
      ! cm1r19:
      j = nj+1
      do while( jflag.lt.0 .and. j.gt.1 )
        j = j-1
        if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
          jflag = j
        endif
      enddo
    ENDIF

  ENDIF  haveit1

#ifdef MPI
      ! check for conflict:
    IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
        (jflag.ge.1.and.jflag.le.nj) )THEN
      IF( iflag.eq.ni .and. pdata(np,prx).eq.xf(iflag+1) .and. nodex.gt.1 .and.  myi.ne.nodex ) iflag = -1
      IF( jflag.eq.nj .and. pdata(np,pry).eq.yf(jflag+1) .and. nodey.gt.1 .and.  myj.ne.nodey ) jflag = -1
    ENDIF
#endif

      myparcel:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                     (jflag.ge.1.and.jflag.le.nj) )THEN

      rkloop:  DO nrkp = 1,2

      IF( nrkp.eq.1 )THEN
        i=iflag
        j=jflag
      ELSE
        iflag = -100
        jflag = -100
        IF(nx.eq.1)THEN
          iflag = 1
        ELSE
          ! cm1r19:
          i = ni+2
          do while( iflag.lt.0 .and. i.gt.0 )
            i = i-1
            if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
              iflag = i
            endif
          enddo
        ENDIF
        IF(axisymm.eq.1.or.ny.eq.1)THEN
          jflag = 1
        ELSE
          do j=0,nj+1
            if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) ) jflag=j
          enddo
          ! cm1r19:
          j = nj+2
          do while( jflag.lt.0 .and. j.gt.0 )
            j = j-1
            if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
              jflag = j
            endif
          enddo
        ENDIF
        i=iflag
        j=jflag
      ENDIF

        IF(debug)THEN
        if( i.lt.0 .or. i.gt.(ni+1) .or. j.lt.0 .or. j.gt.(nj+1) )then
          print *,'  myid,i,j = ',myid,i,j
          print *,'  x,x1     = ',x3d,pdata(np,prx)
          print *,'  y,y1     = ',y3d,pdata(np,pry)
          do i=0,ni+1
            print *,i,abs(xh(i)-x3d),0.5*dx*ruh(i)
          enddo
          do j=0,nj+1
            print *,j,abs(yh(j)-y3d),0.5*dy*rvh(j)
          enddo
          call stopcm1
        endif
        ENDIF

        kflag = 1
        do while( z3d.ge.zf(iflag,jflag,kflag+1) )
          kflag = kflag+1
        enddo

        IF(debug)THEN
        if( kflag.le.0 .or. kflag.ge.(nk+1) )then
          print *,myid,nrkp
          print *,iflag,jflag,kflag
          print *,pdata(np,prx),pdata(np,pry),pdata(np,prz)
          print *,x3d,y3d,z3d
          print *,uval,vval,wval
          print *,zf(iflag,jflag,kflag),z3d,zf(iflag,jflag,kflag+1)
          print *,'  16667 '
          call stopcm1
        endif
        ENDIF

!----------------------------------------------------------------------
!  Data on u points

        i=iflag
        j=jflag
        k=kflag

        if( y3d.lt.yh(j) )then
          j=j-1
        endif
        if( z3d.lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

        rx = ( x3d-xf(i) )/( xf(i+1)-xf(i) )
        ry = ( y3d-yh(j) )/( yh(j+1)-yh(j) )
        rz = ( z3d-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )

        ! saveit:
        rxu = rx
        rys = ry
        rzs = rz

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.0 .or. i.gt.(ni+1)   .or.        &
            j.lt.-1 .or. j.gt.(nj+1)   .or.       &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  13333a: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xf1,x3d,xf2 = ',xf(i),x3d,xf(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),y3d,yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),z3d,zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni+1,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,ua,uval)

!----------------------------------------------------------------------
!  Data on v points

        i=iflag
        j=jflag
        k=kflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( z3d.lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

        rx = ( x3d-xh(i) )/( xh(i+1)-xh(i) )
        ry = ( y3d-yf(j) )/( yf(j+1)-yf(j) )
!!!        rz = ( z3d-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )
        rz = rzs

        ! saveit:
        rxs = rx
        ryv = ry

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.(ni+1)   .or.       &
            j.lt.0 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  23333b: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),x3d,xh(i+1)
          print *,'  yf1,y3d,yh2 = ',yf(j),y3d,yf(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),z3d,zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni,nj+1,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,va,vval)

!----------------------------------------------------------------------
!  uv for parcels below lowest model level:

      IF( bbc.eq.3 )THEN
        ! semi-slip lower boundary condition:
        if( z3d.lt.zh(1,1,1) )then
          ! re-calculate velocities if parcel is below lowest model level:
          !------
          ! u at lowest model level:
          i=iflag
          j=jflag
          if( y3d.lt.yh(j) )then
            j=j-1
          endif
          call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 1, 0, 1, 0,ua(ib,jb,1),uval)
          !------
          ! v at lowest model level:
          i=iflag
          j=jflag
          if( x3d.lt.xh(i) )then
            i=i-1
          endif
          call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 1, 0, 1,va(ib,jb,1),vval)
          !------
          ! z0:
          i=iflag
          j=jflag
          if( x3d.lt.xh(i) )then
            i=i-1
          endif
          if( y3d.lt.yh(j) )then
            j=j-1
          endif
          call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,znt,z0)
          !------
          ! get u,v from (neutral) log-layer equation:
          rznt = 1.0/z0
          var = alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
          uval = uval*var
          vval = vval*var
        endif
      ENDIF

!----------------------------------------------------------------------
!  Data on w points

        i=iflag
        j=jflag
        k=kflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( y3d.lt.yh(j) )then
          j=j-1
        endif

!!!        rx = ( x3d-xh(i) )/( xh(i+1)-xh(i) )
!!!        ry = ( y3d-yh(j) )/( yh(j+1)-yh(j) )
        rx = rxs
        ry = rys
        rz = ( z3d-zf(iflag,jflag,k) )/( zf(iflag,jflag,k+1)-zf(iflag,jflag,k) )

        ! saveit:
        rzw = rz

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.(ni+1)   .or.       &
            j.lt.-1 .or. j.gt.(nj+1)   .or.       &
            k.lt.1 .or. k.gt.nk                   )then
          print *
          print *,'  43333a: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),x3d,xh(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),y3d,yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zf(iflag,jflag,k),z3d,zf(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,wa,wval)

!-----------------------------------------------------
!  Update parcel positions:
!-----------------------------------------------------

      ! RK2 scheme:
      IF(nrkp.eq.1)THEN
        IF(nx.eq.1)THEN
          x3d=0.0
        ELSE
          x3d=pdata(np,prx)+dt*uval
        ENDIF
        IF(axisymm.eq.1.or.ny.eq.1)THEN
          y3d=0.0
        ELSE
          y3d=pdata(np,pry)+dt*vval
        ENDIF
        z3d=pdata(np,prz)+dt*wval
        uu1=uval
        vv1=vval
        ww1=wval
      ELSE
        IF(nx.eq.1)THEN
          x3d=0.0
        ELSE
          x3d=pdata(np,prx)+dt2*(uu1+uval)
        ENDIF
        IF(axisymm.eq.1.or.ny.eq.1)THEN
          y3d=0.0
        ELSE
          y3d=pdata(np,pry)+dt2*(vv1+vval)
        ENDIF
        z3d=pdata(np,prz)+dt2*(ww1+wval)
      ENDIF

        IF( z3d.lt.0.0 )THEN
          print *,'  parcel is below surface:  np,x3d,y3d,z3d = ',np,x3d,y3d,z3d
          z3d=1.0e-6
        ENDIF
        z3d=min(z3d,maxz)

      ENDDO  rkloop

!-----------------------------------------------------
!  Account for boundary conditions (if necessary)
!-----------------------------------------------------

        ! New for cm1r17:  if parcel exits domain,
        ! just assume periodic lateral boundary conditions
        ! (no matter what actual settings are for wbc,ebc,sbc,nbc)

        if(x3d.lt.minx)then
          x3d=x3d+(maxx-minx)
        endif
        if(x3d.gt.maxx)then
          x3d=x3d-(maxx-minx)
        endif

        if( (y3d.gt.maxy).and.(axisymm.ne.1).and.(ny.ne.1) )then
          y3d=y3d-(maxy-miny)
        endif
        if( (y3d.lt.miny).and.(axisymm.ne.1).and.(ny.ne.1) )then
          y3d=y3d+(maxy-miny)
        endif

        pdata(np,prx)=x3d
        pdata(np,pry)=y3d
        pdata(np,prz)=z3d

#ifdef MPI

      ELSE

        ! set to really small number (so we can use the allreduce command below)
        pdata(np,prx) = -1.0e30
        pdata(np,pry) = -1.0e30
        pdata(np,prz) = -1.0e30
#endif

      ENDIF  myparcel

    ENDDO  nploop

!----------------------------------------------------------------------
!  communicate data  (for MPI runs)
#ifdef MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE,pdata(1,1),3*nparcels,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif
!----------------------------------------------------------------------

      end subroutine parcel_driver


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine parcel_interp(dt,xh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,     &
                               zh,mh,rmh,zf,mf,znt,ust,c1,c2,          &
                               pi0,th0,thv0,qv0,qc0,qi0,rth0,          &
                               dum1,dum2,dum3,dum4,zv  ,qt  ,prs,rho,  &
                               dum7,dum8,buoy,vpg  ,                   &
                               u3d,v3d,w3d,pp3d,th   ,t     ,th3d,q3d, &
                               kmh,kmv,khh,khv,tke3d,pt3d,pdata,       &
                               tdiag,qdiag,                            &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,        &
                               nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,reqs_p, &
                               judiag, jvdiag, jwdiag, mtime, rsttime)    ! DAHL 20181204: Added diagnostic arrays and the time variables

      use input
      use constants
      use cm1libs , only : rslf,rsif
      use bc_module
      use comm_module
#ifdef MPI
      use mpi
#endif
      implicit none

!-----------------------------------------------------------------------
!  This subroutine interpolates model information to the parcel locations
!  (diagnostic only ... not used for model integration)
!-----------------------------------------------------------------------

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je) :: znt,ust
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,th0,thv0,qv0,qc0,qi0,rth0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,zv,qt,prs,rho
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum7,dum8
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: w3d,buoy,vpg
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pp3d,th3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: th,t
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tke3d
      real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pt3d
      real, intent(inout), dimension(nparcels,npvals) :: pdata
      real, intent(inout) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
      real, intent(inout) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer, intent(inout), dimension(rmp) :: reqs_p

      integer :: n,np,i,j,k,iflag,jflag,kflag
      real :: tem,tem1
      real :: uval,vval,wval,rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8,wsum
      real :: rxu,ryv,rzw,rxs,rys,rzs
      real :: x3d,y3d,z3d,z0,rznt,var

      logical, parameter :: debug = .false.

! DAHL 20181204

REAL, INTENT(INOUT) :: &
  judiag(ib:ie+1, jb:je  , kb:ke, nudiag), &
  jvdiag(ib:ie  , jb:je+1, kb:ke, nvdiag), &
  jwdiag(ib:ie  , jb:je  , kb:ke+1, nwdiag)

DOUBLE PRECISION, INTENT(IN) :: &
  mtime, rsttime  ! rsttime = const; rsttim updated during the run

! Local variables
! ---------------

! DAHL

REAL :: & 
  dudx(ib:ie,   jb:je,   kb:ke), &
  dudy(ib:ie+1, jb+1:je, kb:ke), &
  dudz(ib:ie+1, jb:je,   kb+1:ke), &

  dvdx(ib+1:ie,jb:je+1,  kb:ke), &
  dvdy(ib:ie,  jb:je,    kb:ke), &
  dvdz(ib:ie,  jb:je+1,  kb+1:ke), &

  dwdx(ib+1:ie,jb:je,   kb:ke+1), &
  dwdy(ib:ie,  jb+1:je, kb:ke+1), &
  dwdz(ib:ie,  jb:je,   kb:ke), &

  xi  (ib:ie,   jb+1:je, kb+1:ke), &
  eta (ib+1:ie, jb:je,   kb+1:ke), &
  zeta(ib+1:ie, jb+1:je, kb:ke)

! Velocity tendencies

REAL :: &
  sgs_u (ib:ie+1,jb:je,   kb:ke),  &
  impl_u(ib:ie+1,jb:je,   kb:ke),  &
  expl_u(ib:ie+1,jb:je,   kb:ke),  &
  pgf_u (ib:ie+1,jb:je,   kb:ke),  &
  cor_u (ib:ie+1,jb:je,   kb:ke),  &

  sgs_v (ib:ie,  jb:je+1, kb:ke),  &
  impl_v(ib:ie,  jb:je+1, kb:ke),  &
  expl_v(ib:ie,  jb:je+1, kb:ke),  &
  pgf_v (ib:ie,  jb:je+1, kb:ke),  &
  cor_v (ib:ie,  jb:je+1, kb:ke),  &

  sgs_w (ib:ie,  jb:je,   kb:ke+1), &
  impl_w(ib:ie,  jb:je,   kb:ke+1), &
  expl_w(ib:ie,  jb:je,   kb:ke+1), &
  pgf_w (ib:ie,  jb:je,   kb:ke+1)

! Turbulent velocity tendency gradients

REAL :: &
  sgs_uy  (ib:ie+1, jb+1:je, kb:ke  ), &
  impl_uy (ib:ie+1, jb+1:je, kb:ke  ), &
  expl_uy (ib:ie+1, jb+1:je, kb:ke  ), &
  pgf_uy  (ib:ie+1, jb+1:je, kb:ke  ), &
  cor_uy  (ib:ie+1, jb+1:je, kb:ke  ), &

  sgs_uz (ib:ie+1, jb:je,   kb+1:ke), &
  impl_uz(ib:ie+1, jb:je,   kb+1:ke), &
  expl_uz(ib:ie+1, jb:je,   kb+1:ke), &
  pgf_uz (ib:ie+1, jb:je,   kb+1:ke), &
  cor_uz (ib:ie+1, jb:je,   kb+1:ke), &

  sgs_vx (ib+1:ie, jb:je+1, kb:ke  ), &
  impl_vx(ib+1:ie, jb:je+1, kb:ke  ), &
  expl_vx(ib+1:ie, jb:je+1, kb:ke  ), &
  pgf_vx (ib+1:ie, jb:je+1, kb:ke  ), &
  cor_vx (ib+1:ie, jb:je+1, kb:ke  ), &

  sgs_vz (ib:ie,   jb:je+1, kb+1:ke), &
  impl_vz(ib:ie,   jb:je+1, kb+1:ke), &
  expl_vz(ib:ie,   jb:je+1, kb+1:ke), &
  pgf_vz (ib:ie,   jb:je+1, kb+1:ke), &
  cor_vz (ib:ie,   jb:je+1, kb+1:ke), &

  sgs_wx  (ib+1:ie, jb:je,   kb:ke+1), &
  impl_wx (ib+1:ie, jb:je,   kb:ke+1), &
  expl_wx(ib+1:ie, jb:je,   kb:ke+1), &
  pgf_wx  (ib+1:ie, jb:je,   kb:ke+1), &

  sgs_wy (ib:ie,   jb+1:je, kb:ke+1), &
  impl_wy(ib:ie,   jb+1:je, kb:ke+1), &
  expl_wy(ib:ie,   jb+1:je, kb:ke+1), &
  pgf_wy (ib:ie,   jb+1:je, kb:ke+1)

REAL ::     &
  psgs_uy,  &  
  pimpl_uy, &
  pexpl_uy, &
  ppgf_uy,  &
  pcor_uy,  &

  psgs_uz,  &
  pimpl_uz, &
  pexpl_uz, &
  ppgf_uz,  &
  pcor_uz,  &

  psgs_vx,  &
  pimpl_vx, &
  pexpl_vx, &
  ppgf_vx,  &
  pcor_vx,  &

  psgs_vz,  &
  pimpl_vz, &
  pexpl_vz, &
  ppgf_vz,  &
  pcor_vz,  &

  psgs_wx,  &
  pimpl_wx, &
  pexpl_wx, &
  ppgf_wx,  &

  psgs_wy,  &
  pimpl_wy, &
  pexpl_wy, &
  ppgf_wy 


INTEGER :: &
  istat, xstag, ystag, zstag

REAL :: &
  x_xi(ib:ie),   &
  y_xi(jb+1:je), &
  z_xi(kb+1:ke), &

  x_eta(ib+1:ie), &
  y_eta(jb:je)  , &
  z_eta(kb+1:ke), &

  x_zeta(ib+1:ie), &
  y_zeta(jb+1:je), &
  z_zeta(kb:ke)

REAL     ::      &
  x_dudx(ib:ie), &
  y_dudx(jb:je), &
  z_dudx(kb:ke), &
  x_dudy(ib:ie+1), &
  y_dudy(jb+1:je), &
  z_dudy(kb:ke),   &

  x_dudz(ib:ie+1), &
  y_dudz(jb:je),   &
  z_dudz(kb+1:ke), &

  x_dvdx(ib+1:ie), &
  y_dvdx(jb:je+1), &
  z_dvdx (kb:ke),  &

  x_dvdy(ib:ie), &
  y_dvdy(jb:je), &
  z_dvdy(kb:ke), &

  x_dvdz(ib:ie),   &
  y_dvdz(jb:je+1), &
  z_dvdz(kb+1:ke), &

  x_dwdx(ib+1:ie), &
  y_dwdx(jb:je),   &
  z_dwdx(kb:ke+1), &

  x_dwdy(ib:ie),   &
  y_dwdy(jb+1:je), &
  z_dwdy(kb:ke+1),   &

  x_dwdz(ib:ie), &
  y_dwdz(jb:je), &
  z_dwdz(kb:ke)

REAL :: &
  x_dphidx(ib+1:ie), &
  y_dphidx(jb  :je), &
  z_dphidx(kb  :ke), &

  x_dphidy(ib  :ie), &
  y_dphidy(jb+1:je), &
  z_dphidy(kb  :ke), &

  x_dphidz(ib  :ie), &
  y_dphidz(jb  :je), &
  z_dphidz(kb+1:ke), &

  x_d2phixx(ib+1:ie-1), &
  y_d2phixx(jb  :je  ), &
  z_d2phixx(kb  :ke  ), &

  x_d2phixy(ib+1:ie), &
  y_d2phixy(jb+1:je), &
  z_d2phixy(kb  :ke), &

  x_d2phiyy(ib  :ie  ), &
  y_d2phiyy(jb+1:je-1), &
  z_d2phiyy(kb  : ke ), &

  x_d2phixz(ib+1:ie), &
  y_d2phixz(jb  :je), &
  z_d2phixz(kb+1:ke), &
  x_d2phizz(ib  :ie  ), &
  y_d2phizz(jb  :je  ), &
  z_d2phizz(kb  :ke  )

REAL :: &

! Exner

  dpidx(ib+1:ie, jb  :je, kb  :ke), &
  dpidy(ib  :ie, jb+1:je, kb  :ke), &
  dpidz(ib  :ie, jb  :je, kb+1:ke), &

  d2pixx(ib+1:ie-1, jb  :je,   kb  :ke), &
  d2pixy(ib+1:ie,   jb+1:je,   kb  :ke), &
  d2piyy(ib  :ie,   jb+1:je-1, kb  :ke), &
  d2pixz(ib+1:ie,   jb:je,     kb+1:ke), &
  d2pizz(ib  :ie,   jb:je,     kb  :ke), &

! Buoyancy

  dbdx(ib+1:ie, jb  :je, kb  :ke), &
  dbdy(ib  :ie, jb+1:je, kb  :ke), &
  dbdz(ib  :ie, jb  :je, kb+1:ke)

INTEGER :: &
  istart, iend, &
  jstart, jend, &
  kstart, kend, &
  ii, jj, kk

! Variables interpolated to parcels

REAL :: &
  tempvar, &

  parcelrho, & ! Density interpolated to parcel location

! Velocity gradients

  pdudx, &
  pdudy, &
  pdudz, &

  pdvdx, &
  pdvdy, &
  pdvdz, &

  pdwdx, &
  pdwdy, &
  pdwdz, &

! Vorticity forcing along trajectories

  pxiten_sgs,   &     
  pxiten_idiff, & 
  pxiten_ediff, &
  pxiten_pgf,   &
  pxiten_cor,   &

  petaten_sgs,   &
  petaten_idiff, &
  petaten_ediff, &
  petaten_pgf,   &
  petaten_cor,   &

  pzetaten_sgs,   &
  pzetaten_idiff, &
  pzetaten_ediff, &
  pzetaten_pgf,   &
  pzetaten_cor,   &

! Velocity variables

  vmagn,         &
  ufrac,         &
  vfrac,         &
  wfrac,         &
  pdvmdx,        &
  pdvmdy,        &
  pdvmdz,        &

  pxi,           &
  peta,          &
  pzeta,         &

  eta_shr,       &
  eta_crv,       &
  eta_conv1,     &
  eta_conv2,     &
  eta_conv3,     &
  eta_crv_bgen1, &
  eta_crv_bgen2, &
  eta_crv_bgen3, &
  eta_shr_bgen,  &

  zeta_shr,      &
  zeta_crv,      &
  zeta_conv1,    &
  zeta_conv2,    &
  omega_cw,      &
  omega_sw,      &
  zeta_shr_tilt, &
  zeta_crv_tilt , &
  zeta_shr_stretch, &
  zeta_crv_stretch, &

  ppsi,            &  ! azimuthal angle

! Turbulent tendencies

  pdutfdy, &
  pdutfdz, &
  pdvtfdx, &
  pdvtfdz, &
  pdwtfdx, &
  pdwtfdy, &

  xitforcing, &
  etatforcing, &
  zetatforcing

REAL :: bb = 1.E-6  ! Small number

INTEGER :: &
  ilbc = 1  ! 1: New/zero flux gradient velicity profile; 2: Traditional zero-gradient profile

! Debug / temporary

REAL :: &
  t1, t2, t3, tick

INTEGER :: kkk, aa

!----------------------------------------------------------------------
!  Get derived variables:

    IF(imoist.eq.1)THEN
      ! with moisture:

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n)
    do k=1,nk

      do j=1,nj
      do i=1,ni
        qt(i,j,k)=q3d(i,j,k,nqv)
      enddo
      enddo
      do n=nql1,nql2
        do j=1,nj
        do i=1,ni
          qt(i,j,k)=qt(i,j,k)+q3d(i,j,k,n)
        enddo
        enddo
      enddo
      IF(iice.eq.1)THEN
        do n=nqs1,nqs2
        do j=1,nj
        do i=1,ni
          qt(i,j,k)=qt(i,j,k)+q3d(i,j,k,n)
        enddo
        enddo
        enddo
      ENDIF
      IF( prth.ge.1 .or. prt.ge.1 .or. prqsl.ge.1 .or. prqsi.ge.1 .or.  prvpg.ge.1 )THEN
        do j=1,nj
        do i=1,ni
          th(i,j,k) = (th0(i,j,k)+th3d(i,j,k))
          t(i,j,k) = th(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))
        enddo
        enddo
      ENDIF
      IF( prb.ge.1 .or. prvpg.ge.1 )THEN
        do j=1,nj
        do i=1,ni
          dum7(i,j,k) = g*( th3d(i,j,k)*rth0(i,j,k)             &
                           +repsm1*(q3d(i,j,k,nqv)-qv0(i,j,k))  &
                           -(qt(i,j,k)-q3d(i,j,k,nqv)-qc0(i,j,k)-qi0(i,j,k))   )
        enddo
        enddo
      ENDIF
      IF( prvpg.ge.1 )THEN
        do j=1,nj
        do i=1,ni
          dum8(i,j,k) = th(i,j,k)*(1.0+reps*q3d(i,j,k,nqv))/(1.0+qt(i,j,k))
        enddo
        enddo
      ENDIF

    enddo

    ELSE
      ! dry:

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
    do k=1,nk

      IF( prth.ge.1 .or. prt.ge.1 .or. prvpg.ge.1 )THEN
        do j=1,nj
        do i=1,ni
          th(i,j,k)= (th0(i,j,k)+th3d(i,j,k))
          t(i,j,k) = th(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))
        enddo
        enddo
      ENDIF
      IF( prb.ge.1 .or. prvpg.ge.1 )THEN
        do j=1,nj
        do i=1,ni
          dum7(i,j,k) = g*( th3d(i,j,k)*rth0(i,j,k) )
        enddo
        enddo
      ENDIF
      IF( prvpg.ge.1 )THEN
        do j=1,nj
        do i=1,ni
          dum8(i,j,k) = th(i,j,k)
        enddo
        enddo
      ENDIF

    enddo

    ENDIF


    IF( prb.ge.1 .or. prvpg.ge.1 )THEN
      do k=2,nk
      do j=1,nj
      do i=1,ni
        buoy(i,j,k) = (c1(1,1,k)*dum7(i,j,k-1)+c2(1,1,k)*dum7(i,j,k))
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_parcels=time_parcels+mytime()
      call prepcornert(buoy,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                            pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      do j=0,nj+1
      do i=0,ni+1
        buoy(i,j,1) = buoy(i,j,2)+(buoy(i,j,3)-buoy(i,j,2))  &
                                 *(  zf(i,j,1)-  zf(i,j,2))  &
                                 /(  zf(i,j,3)-  zf(i,j,2))
        buoy(i,j,nk+1) = buoy(i,j,nk)+(buoy(i,j,nk  )-buoy(i,j,nk-1))  &
                                     *(  zf(i,j,nk+1)-  zf(i,j,nk  ))  &
                                     /(  zf(i,j,nk  )-  zf(i,j,nk-1))
      enddo
      enddo
    ENDIF
    IF( prvpg.ge.1 )THEN
      tem1 = rdz*cp
      ! assuming no terrain:
      do k=2,nk
      do j=1,nj
      do i=1,ni
        vpg(i,j,k) = -tem1*(pp3d(i,j,k)-pp3d(i,j,k-1))*mf(1,1,k)  &
                          *(c2(1,1,k)*dum8(i,j,k)+c1(1,1,k)*dum8(i,j,k-1))
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_parcels=time_parcels+mytime()
      call prepcornert(vpg,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                           pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      ! cmr18:  at top/bottom boundaries, vpg + buoy = 0
      do j=0,nj+1
      do i=0,ni+1
        vpg(i,j,1) = -buoy(i,j,1)
        vpg(i,j,nk+1) = -buoy(i,j,nk+1)
      enddo
      enddo
    ENDIF

    if(timestats.ge.1) time_parcels=time_parcels+mytime()

!----------------------------------------------------------------------
!  get corner info for MPI runs
!  (may not parallelize correctly if this is not done)

#ifdef MPI
      call getcorneru(u3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcu2(u3d)
      call getcornerv(v3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcv2(v3d)
      call getcornerw(w3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcw2(w3d)
#endif

!----------------------------------------------------------------------
!  apply bottom/top boundary conditions:
!  [Note:  for u,v,s the array index (i,j,0) means the surface, ie z=0]
!     (for the parcel subroutines only!)

!$omp parallel do default(shared)  &
!$omp private(i,j)

IF(bbc.eq.1)THEN

! free slip ... extrapolate:

  DO j=jb,je
    do i=ib,ie+1
      u3d     (i,j,0) = cgs1*u3d   (i,j,1)  +cgs2*u3d   (i,j,2)+cgs3*u3d     (i,j,3)
      judiag(i,j,0,:) = cgs1*judiag(i,j,1,:)+cgs2*judiag(i,j,2,:)+cgs3*judiag(i,j,3,:)
    enddo
  ENDDO

  DO j = jb, je+1
    do i=ib,ie
      v3d    (i,j,0)  = cgs1*v3d   (i,j,1)  +cgs2*v3d   (i,j,2)+cgs3*v3d     (i,j,3)
      jvdiag(i,j,0,:) = cgs1*jvdiag(i,j,1,:)+cgs2*jvdiag(i,j,2,:)+cgs3*jvdiag(i,j,3,:)
    enddo
  ENDDO
ELSEIF(bbc.eq.2)THEN

! The LBC for the diagnostic fields are really just place holders.  It's safest
! to discard all parcels that descend below the lowest scalar model level.

      ! no slip:
  DO j = jb, je
    do i = ib, ie+1
        u3d   (i,j,0)   = 0.0
        judiag(i,j,0,:) = judiag(i,j,1,:)
     enddo
  ENDDO

  DO j = jb, je+1
    do i=ib,ie
      v3d(i,j,0) = 0.0
      jvdiag(i,j,0,:) = jvdiag(i,j,1,:)
    enddo
  ENDDO
ELSEIF(bbc.eq.3)THEN
      ! u,v near sfc are determined below using log-layer equations

! Setting the velocity accelerations to constant values below the lowest
! scalar level.  This is a somewhat random choice and it is best to discard
! any parcel that even for some time descends below the lowest scalar level.

  DO j = jb, je
    DO i = ib, ie+1
      judiag(i,j,0,:) =  judiag(i,j,1,:)
     ENDDO
  ENDDO

  DO j = jb, je+1
    DO i = ib, ie
      jvdiag(i,j,0,:) =  jvdiag(i,j,1,:)
    ENDDO
  ENDDO

ENDIF

!----------

  DO j = jb, je+1

    IF(tbc.eq.1)THEN
      ! free slip ... extrapolate:
      IF(j.le.je)THEN
      do i=ib,ie+1
        u3d(i,j,nk+1) = cgt1*u3d(i,j,nk)+cgt2*u3d(i,j,nk-1)+cgt3*u3d(i,j,nk-2)
      enddo
      ENDIF
      do i=ib,ie
        v3d(i,j,nk+1) = cgt1*v3d(i,j,nk)+cgt2*v3d(i,j,nk-1)+cgt3*v3d(i,j,nk-2)
      enddo
    ELSEIF(tbc.eq.2)THEN
      ! no slip:
      IF(j.le.je)THEN
      do i=ib,ie+1
        u3d(i,j,nk+1) = 0.0
      enddo
      ENDIF
      do i=ib,ie
        v3d(i,j,nk+1) = 0.0
      enddo
    ENDIF

!----------

      ! assuming no terrain:
      IF(j.le.je)THEN
      do i=ib,ie
        w3d(i,j,1)    = 0.0
        w3d(i,j,nk+1) = 0.0
      enddo
      ENDIF

  ENDDO

      if(timestats.ge.1) time_parcels=time_parcels+mytime()

      if( prth.ge.1 )then
        call prepcorners(th ,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif
      if( prt.ge.1 )then
        call prepcorners(t  ,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif
      if( prprs.ge.1 )then
        call prepcorners(prs,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif
!      if( prrho.ge.1 )then   ! DAHL always do
        call prepcorners(rho,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
!      endif
      if(prpt1.ge.1)then
        do n=1,npt
          call prepcorners(pt3d(ib,jb,kb,n),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                                            pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
        enddo
      endif
      if( prqv.ge.1 )then
        call prepcorners(q3d(ib,jb,kb,nqv),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                                           pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
      endif
      if( prq1.ge.1 .or. prnc1.ge.1 )then
        do n = 1,numq
          call prepcorners(q3d(ib,jb,kb,n),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                                           pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
        enddo
      endif
      if( prkm.ge.1 )then
        call prepcornert(kmh,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
        call prepcornert(kmv,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
      endif
      if( prkh.ge.1 )then
        call prepcornert(khh,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
        call prepcornert(khv,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
      endif
      if( prtke.ge.1 )then
        call prepcornert(tke3d,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
      endif


    ! DAHL 6/9/2023 Removed this code because it's not needed for the vorticity
    ! diagnostics, and since qd_dbz tends to be zero sometimes, this leads to
    ! segmentation violations.

    !  if( prdbz.ge.1 )then
    !    call prepcorners(qdiag(ibdq,jbdq,kbdq,qd_dbz),  &
    !                         nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
    !                         pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
    !  endif

!----------------------------------------------------------------------

    IF( prqsl.ge.1 )THEN
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum1(i,j,k) = rslf( prs(i,j,k) , t(i,j,k) )
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_parcels=time_parcels+mytime()
      call prepcorners(dum1,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                            pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
    ENDIF
    IF( prqsi.ge.1 )THEN
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum2(i,j,k) = rsif( prs(i,j,k) , t(i,j,k) )
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_parcels=time_parcels+mytime()
      call prepcorners(dum2,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                            pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
    ENDIF

!----------------------------------------------------------------------
!  Get zvort at appropriate C-grid location:
!  (assuming no terrain)
!  cm1r18:  below lowest model level:
!           Use extrapolated velocities for bbc=1,2
!           Use log-layer equations for bbc=3 (see below)

    IF( przv.ge.1)THEN

      do k=0,nk+1
      do j=1,nj+1
      do i=1,ni+1
        zv(i,j,k) = (v3d(i,j,k)-v3d(i-1,j,k))*rdx*uf(i)   &
                   -(u3d(i,j,k)-u3d(i,j-1,k))*rdy*vf(j)
      enddo
      enddo
      enddo

    ENDIF

! DAHL 4/4/18

!zzh        = zh
!zzh(:,:,0) = -zzh(:,:,1) 

!IF (myid == 1) THEN
!  DO k = kb, ke
! print *, 'AA', k, zzh(1,1,k)
!ENDDO
!ENDIF
! Adding forcing functions (new arrays combining vertical and
! horizontal diffusion/mixig tendencies)

sgs_u = 0.0
sgs_v = 0.0
sgs_w = 0.0

impl_u = 0.0
impl_v = 0.0
impl_w = 0.0

expl_u = 0.0
expl_v = 0.0
expl_w = 0.0

pgf_u = 0.0
pgf_v = 0.0
pgf_w = 0.0

cor_u = 0.0
cor_v = 0.0

sgs_u(:,:,:) = judiag(:,:,:, ud_hturb) + judiag(:,:,:, ud_vturb)  
sgs_v (:,:,:)=jvdiag(:,:,:, vd_hturb)  + jvdiag(:,:,:, vd_vturb) 
sgs_w(:,:,:) = jwdiag(:,:,:, wd_hturb) + jwdiag(:,:,:, wd_vturb)

! A little error management

   IF (hadvordrv .NE. vadvordrv) THEN
  WRITE(*,*) &
  'For the vorticity diagnostic along he trajectories, the order of the horizontal advection scheme must match order of the  vertical advection scheme'
  WRITE(*,*) 'Stopping CM1'
  CALL stopcm1
ENDIF

IF (vadvordrv.eq.3 .or. vadvordrv.eq.5 .or. vadvordrv.eq.7 .or. vadvordrv.eq.9 .or. advwenov.ge.1) THEN

  impl_u(:,:,:) = judiag(:,:,:, ud_hidiff) + judiag(:,:,:, ud_vidiff) 
  impl_v(:,:,:) =  jvdiag(:,:,:, vd_hidiff) + jvdiag(:,:,:, vd_vidiff)
  impl_w(:,:,:) = jwdiag(:,:,:, wd_hidiff) + jwdiag(:,:,:, wd_vidiff)

ENDIF

IF (idiff >=1) THEN
   expl_u(:,:,:) = judiag(:,:,:, ud_hediff)
   expl_v(:,:,:) =jvdiag(:,:,:, vd_hediff)
   expl_w(:,:,:) =jwdiag(:,:,:, wd_hediff)
  if( difforder.eq.2 )then
    expl_u(:,:,:) = judiag(:,:,:, ud_hediff) + judiag(:,:,:, ud_vediff)
    expl_v(:,:,:) =jvdiag(:,:,:, vd_hediff) + jvdiag(:,:,:, vd_vediff) 
    expl_w(:,:,:) =jwdiag(:,:,:, wd_hediff) + jwdiag(:,:,:, wd_vediff)
  endif
ENDIF
              
pgf_u(:,:,:) =judiag(:,:,:, ud_pgrad)
pgf_v(:,:,:) =jvdiag(:,:,:, vd_pgrad)
pgf_w(:,:,:) =jwdiag(:,:,:, wd_pgrad) + jwdiag(:,:,:, wd_buoy)
                 
IF (icor == 1) THEN
  cor_u(:,:,:) = judiag(:,:,:, ud_cor)
  cor_v(:,:,:) = jvdiag(:,:,:, vd_cor)                          
ENDIF

! Since the interpolator is called before the solver AT the restart time,
! we need to initialize the forcing functions with zero for the first time step.
! Re-initializing seems more efficient than having the above within an
! IF-branch.

IF (mtime < (rsttime + bb) .AND. mtime > (rsttime - bb)) THEN

  sgs_u = 0.0
  sgs_v = 0.0
  sgs_w = 0.0

  impl_u = 0.0
  impl_v = 0.0
  impl_w = 0.0

  expl_u = 0.0
  expl_v = 0.0
  expl_w = 0.0

  pgf_u = 0.0
  pgf_v = 0.0
  pgf_w = 0.0

  cor_u = 0.0
  cor_v = 0.0

ENDIF

!IF (myid == 44) THEN
!  print *, 'DD1', MAXVAL(sgs_u), MAXVAL(impl_u), MAXVAL(expl_u), MAXVAL(pgf_u), MAXVAL(cor_u)
!  print *, mtime, rsttime
!ENDIF

! zh is the original zh, which is zero for k = 0.

CALL delv       (ilbc, rmh, mf, xh, yh, zh, zf,                         &

                 ! Velocity gradients, velocity components, vorticity components

                 dudx, dudy, dudz, dvdx, dvdy, dvdz,                     &
                 dwdx, dwdy, dwdz, u3d, v3d, w3d,  xi, eta, zeta,        &

                 ! 3D accelerations

                 sgs_u,  sgs_v,  sgs_w,                                  &
                 impl_u, impl_v, impl_w,                                 &
                 expl_u, expl_v, expl_w,                                 &
                 pgf_u,  pgf_v,  pgf_w,                                  &
                 cor_u,  cor_v,                                          &

                 ! Antisymmetric derivatives of 3D accelerations

                 sgs_uy, impl_uy, expl_uy, pgf_uy, cor_uy,               &
                 sgs_uz, impl_uz, expl_uz, pgf_uz, cor_uz,               &
                 sgs_vx, impl_vx, expl_vx, pgf_vx, cor_vx,               &
                 sgs_vz, impl_vz, expl_vz, pgf_vz, cor_vz,               &
                 sgs_wx, impl_wx, expl_wx, pgf_wx,                       &
                 sgs_wy, impl_wy, expl_wy, pgf_wy,                       &

                 ! Grids associated with the various gradients

                 x_dudx, y_dudx, z_dudx, x_dudy, y_dudy, z_dudy,         &
                 x_dudz, y_dudz, z_dudz, x_dvdx, y_dvdx, z_dvdx,         &
                 x_dvdy, y_dvdy, z_dvdy, x_dvdz, y_dvdz, z_dvdz,         &
                 x_dwdx, y_dwdx, z_dwdx, x_dwdy, y_dwdy, z_dwdy,         &
                 x_dwdz, y_dwdz, z_dwdz, x_xi, y_xi, z_xi,               &
                 x_eta, y_eta, z_eta, x_zeta, y_zeta, z_zeta)



!IF (myid == 26) THEN
!  DO k = kb, ke
!     print *, 'LL1', k, sgs_u(72, 40, k), zh(3,3,k)
!  ENDDO
! print *, ''

!  DO k = kb+1, ke
!     print *, 'LL2', k, sgs_uz(72, 40, k), z_dudz(k)
!  ENDDO
!ENDIF

! rho = prs / (287.04 * th * (pi0 + ppi)) ! I think this is no longer needed

!----------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!----------------------------------------------------------------------
!  Loop through all parcels:  if you have it, get interpolated info:

    nploop2:  &
    DO np=1,nparcels

      x3d = pdata(np,prx)
      y3d = pdata(np,pry)
      z3d = pdata(np,prz)

      iflag = -100
      jflag = -100
      kflag = 0

  ! cm1r19:  skip if we already know this processor doesnt have this parcel
  haveit2:  &
  IF( x3d.ge.xf(1) .and. x3d.le.xf(ni+1) .and.  &
      y3d.ge.yf(1) .and. y3d.le.yf(nj+1) )THEN

    IF(nx.eq.1)THEN
      iflag = 1
    ELSE
      ! cm1r19:
      i = ni+1
      do while( iflag.lt.0 .and. i.gt.1 )
        i = i-1
        if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
          iflag = i
        endif
      enddo
    ENDIF

    IF(axisymm.eq.1.or.ny.eq.1)THEN
      jflag = 1
    ELSE
      ! cm1r19:
      j = nj+1
      do while( jflag.lt.0 .and. j.gt.1 )
        j = j-1
        if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
          jflag = j
        endif
      enddo
    ENDIF

  ENDIF  haveit2

#ifdef MPI
      ! check for conflict:
    IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
        (jflag.ge.1.and.jflag.le.nj) )THEN
      IF( iflag.eq.ni .and. pdata(np,prx).eq.xf(iflag+1) .and. nodex.gt.1 .and.  myi.ne.nodex ) iflag = -1
      IF( jflag.eq.nj .and. pdata(np,pry).eq.yf(jflag+1) .and. nodey.gt.1 .and.  myj.ne.nodey ) jflag = -1
    ENDIF
#endif

      myprcl:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                   (jflag.ge.1.and.jflag.le.nj) )THEN

        i=iflag
        j=jflag

        kflag = 1
        do while( pdata(np,prz).ge.zf(iflag,jflag,kflag+1) )
          kflag = kflag+1
        enddo

        x3d = pdata(np,prx)
        y3d = pdata(np,pry)
        z3d = pdata(np,prz)

IF (myid == 1) print *, myid, 'LLL', x3d, y3d, z3d
!=============================================================================
! DAHL: Adding 21 new interpolators for 21 different grids d/dx_j(u_i)
!       and omega.
!
!                    *** NOTE ***
!
! If the parcels reach the lateral edges of the domain, the interpolator may
! encounter segmentation errors.  This can be fixed by assigning boundary
! values to the fields (not all are affected; only those fields that got
! smaller in size after the differentiation).
!
!=============================================================================



!-----------------------------------------------------------------------------
! 1. dudx
!-----------------------------------------------------------------------------

        istart = ib
        iend   = ie

        jstart = jb
        jend   = je

        kstart = kb
        kend   = ke

        ! Not fast, but safe

        DO ii = istart, iend-1
          IF ( (pdata(np,1) >= x_dudx(ii)) .AND. (pdata(np,1) < x_dudx(ii+1)) ) THEN
            i = ii
            EXIT
          ENDIF
        ENDDO

        DO jj = jstart, jend-1
          IF ( (pdata(np,2) >= y_dudx(jj)) .AND. (pdata(np,2) < y_dudx(jj+1)) ) THEN
            j = jj
            EXIT
          ENDIF
        ENDDO

        DO kk = kstart, kend-1
          IF ( (pdata(np,3) >= z_dudx(kk)) .AND. (pdata(np,3) < z_dudx(kk+1)) ) THEN
            k = kk
            EXIT
          ENDIF
        ENDDO

        rx = ( pdata(np,1)-x_dudx(i)) * rdx
        ry = ( pdata(np,2)-y_dudx(j)) * rdy
        rz = ( pdata(np,3)-z_dudx(k)) / (z_dudx(k+1) - z_dudx(k))

        ! Sanity check

       IF (rx > 1.0 .OR. rx < 0.0 .OR. &
            ry > 1.0 .OR. ry < 0.0 .OR. &
            rz > 1.0 .OR. rz < 0.0 ) THEN

           WRITE(*,*) '    *** Interpolation error: dudx ***'
           WRITE(*,*) 'rx, ry, rz: ', rx, ry, rz
           WRITE(*,*) 'parcel: ', np
           WRITE(*,*) 'Aborting ...'
           CALL stopcm1

        ENDIF

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dudx,pdudx)

!-----------------------------------------------------------------------------
! 2. dudy
!-----------------------------------------------------------------------------

        istart = ib
        iend   = ie+1

        jstart = jb+1
        jend   = je

        kstart = kb
        kend   = ke

        ! Not fast, but safe

        DO ii = istart, iend-1
          IF ( (pdata(np,1) >= x_dudy(ii)) .AND. (pdata(np,1) < x_dudy(ii+1)) ) THEN
            i = ii
            EXIT
          ENDIF
        ENDDO

        DO jj = jstart, jend-1
          IF ( (pdata(np,2) >= y_dudy(jj)) .AND. (pdata(np,2) < y_dudy(jj+1)) ) THEN
            j = jj
            
            EXIT
          ENDIF
        ENDDO

        DO kk = kstart, kend-1
          IF ( (pdata(np,3) >= z_dudy(kk)) .AND. (pdata(np,3) < z_dudy(kk+1)) ) THEN
            k = kk
            EXIT
          ENDIF
        ENDDO

        rx = ( pdata(np,1)-x_dudy(i)) * rdx
        ry = ( pdata(np,2)-y_dudy(j)) * rdy
        rz = ( pdata(np,3)-z_dudy(k)) / (z_dudy(k+1) - z_dudy(k))

        ! Sanity check

        IF (rx > 1.0 .OR. rx < 0.0 .OR. &
            ry > 1.0 .OR. ry < 0.0 .OR. &
            rz > 1.0 .OR. rz < 0.0 ) THEN

           WRITE(*,*) '    *** Interpolation error: dudy ***'
           WRITE(*,*) 'rx, ry, rz: ', rx, ry, rz
           WRITE(*,*) 'parcel: ', np
           WRITE(*,*) 'Aborting ...'
           CALL stopcm1

        ENDIF
        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dudy,pdudy)

        ! u-tendencies

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,sgs_uy, psgs_uy)


 !       print *, myid, 'LL',  np, psgs_uy , pdata(np,2), j  


        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,impl_uy, pimpl_uy)


        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,expl_uy, pexpl_uy)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,pgf_uy, ppgf_uy)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,cor_uy, pcor_uy)

!-----------------------------------------------------------------------------
! 3. dudz
!-----------------------------------------------------------------------------

        istart = ib
        iend   = ie+1

        jstart = jb
        jend   = je

        kstart = kb+1
        kend   = ke

        ! Not fast, but safe

        DO ii = istart, iend-1
          IF ( (pdata(np,1) >= x_dudz(ii)) .AND. (pdata(np,1) < x_dudz(ii+1)) ) THEN
            i = ii
            EXIT
          ENDIF
        ENDDO

        DO jj = jstart, jend-1
          IF ( (pdata(np,2) >= y_dudz(jj)) .AND. (pdata(np,2) < y_dudz(jj+1)) ) THEN
            j = jj
            EXIT
          ENDIF
        ENDDO

        DO kk = kstart, kend-1
          IF ( (pdata(np,3) >= z_dudz(kk)) .AND. (pdata(np,3) < z_dudz(kk+1)) ) THEN
            k = kk
            EXIT
          ENDIF
        ENDDO

        rx = ( pdata(np,1)-x_dudz(i)) * rdx
        ry = ( pdata(np,2)-y_dudz(j)) * rdy
        rz = ( pdata(np,3)-z_dudz(k)) / (z_dudz(k+1) - z_dudz(k))

        ! Sanity check

        IF (rx > 1.0 .OR. rx < 0.0 .OR. &
            ry > 1.0 .OR. ry < 0.0 .OR. &
            rz > 1.0 .OR. rz < 0.0 ) THEN

           WRITE(*,*) '    *** Interpolation error: dudz ***'
           WRITE(*,*) 'rx, ry, rz: ', rx, ry, rz
           WRITE(*,*) 'parcel: ', np
           WRITE(*,*) 'Aborting ...'
           CALL stopcm1

        ENDIF

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dudz,pdudz)

        ! d/dz: u-momentum forcing

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,sgs_uz,psgs_uz)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,impl_uz,pimpl_uz)


        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,expl_uz, pexpl_uz)


        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,pgf_uz,ppgf_uz)


        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,cor_uz,pcor_uz)

!-----------------------------------------------------------------------------
! 4. dvdx
!-----------------------------------------------------------------------------

!       DO aa = 1, 40 
!         pdata(np,1) = 20500.00 + 100.0 * aa


        istart = ib+1
        iend   = ie

        jstart = jb
        jend   = je+1

        kstart = kb
        kend   = ke

        DO ii = istart, iend-1
          IF ( (pdata(np,1) >= x_dvdx(ii)) .AND. (pdata(np,1) < x_dvdx(ii+1))) THEN
            i = ii
            EXIT
          ENDIF
        ENDDO

        DO jj = jstart, jend-1
          IF ( (pdata(np,2) >= y_dvdx(jj)) .AND. (pdata(np,2) < y_dvdx(jj+1)) ) THEN
            j = jj
            EXIT
          ENDIF
        ENDDO

        DO kk = kstart, kend-1
          IF ( (pdata(np,3) >= z_dvdx(kk)) .AND. (pdata(np,3) < z_dvdx(kk+1)) ) THEN
            k = kk
            EXIT
          ENDIF
        ENDDO
        rx = ( pdata(np,1)-x_dvdx(i)) * rdx
        ry = ( pdata(np,2)-y_dvdx(j)) * rdy
        rz = ( pdata(np,3)-z_dvdx(k)) / (z_dvdx(k+1) - z_dvdx(k))

        ! Sanity check

        IF (rx > 1.0 .OR. rx < 0.0 .OR. &
            ry > 1.0 .OR. ry < 0.0 .OR. &
            rz > 1.0 .OR. rz < 0.0 ) THEN

           WRITE(*,*) '    *** Interpolation error: dvdx  ***'
           WRITE(*,*) 'rx, ry, rz: ', rx, ry, rz
           WRITE(*,*) 'parcel: ', np
           WRITE(*,*) 'Aborting ...'
           CALL stopcm1

        ENDIF

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dvdx,pdvdx)

        ! x-gradients of Momentum forcing

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,sgs_vx,psgs_vx)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,impl_vx,pimpl_vx)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,expl_vx,pexpl_vx)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,pgf_vx,ppgf_vx)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,cor_vx,pcor_vx)

!-----------------------------------------------------------------------------
! 5. dvdy
!-----------------------------------------------------------------------------

        istart = ib
        iend   = ie

        jstart = jb
        jend   = je

        kstart = kb
        kend   = ke

        DO ii = istart, iend-1
          IF ( (pdata(np,1) >= x_dvdy(ii)) .AND. (pdata(np,1) < x_dvdy(ii+1)) ) THEN
            i = ii
            EXIT
          ENDIF
        ENDDO

        DO jj = jstart, jend-1
          IF ( (pdata(np,2) >= y_dvdy(jj)) .AND. (pdata(np,2) < y_dvdy(jj+1)) ) THEN
            j = jj
            EXIT
          ENDIF
        ENDDO

        DO kk = kstart, kend-1
          IF ( (pdata(np,3) >= z_dvdy(kk)) .AND. (pdata(np,3) < z_dvdy(kk+1)) ) THEN
            k = kk
            EXIT
          ENDIF
        ENDDO

        rx = ( pdata(np,1)-x_dvdy(i)) * rdx
        ry = ( pdata(np,2)-y_dvdy(j)) * rdy
        rz = ( pdata(np,3)-z_dvdy(k)) / (z_dvdy(k+1) - z_dvdy(k))

        ! Sanity check

        IF (rx > 1.0 .OR. rx < 0.0 .OR. &
            ry > 1.0 .OR. ry < 0.0 .OR. &
            rz > 1.0 .OR. rz < 0.0 ) THEN

           WRITE(*,*) '    *** Interpolation error: dvdy ***'
           WRITE(*,*) 'rx, ry, rz: ', rx, ry, rz
           WRITE(*,*) 'parcel: ', np
           WRITE(*,*) 'Aborting ...'
           CALL stopcm1

        ENDIF

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dvdy,pdvdy)

!-----------------------------------------------------------------------------
! 6. dvdz
!-----------------------------------------------------------------------------

        istart = ib
        iend   = ie

        jstart = jb
        jend   = je+1

        kstart = kb+1
        kend   = ke

        DO ii = istart, iend-1
          IF ( (pdata(np,1) >= x_dvdz(ii)) .AND. (pdata(np,1) < x_dvdz(ii+1))) THEN
            i = ii
            EXIT
          ENDIF
        ENDDO

        DO jj = jstart, jend-1
          IF ( (pdata(np,2) >= y_dvdz(jj)) .AND. (pdata(np,2) < y_dvdz(jj+1)) ) THEN
            j = jj
            EXIT
          ENDIF
        ENDDO

        DO kk = kstart, kend-1
          IF ( (pdata(np,3) >= z_dvdz(kk)) .AND. (pdata(np,3) < z_dvdz(kk+1)) ) THEN
            k = kk
            EXIT
          ENDIF
        ENDDO
        rx = ( pdata(np,1)-x_dvdz(i)) * rdx
        ry = ( pdata(np,2)-y_dvdz(j)) * rdy
        rz = ( pdata(np,3)-z_dvdz(k)) / (z_dvdz(k+1) - z_dvdz(k))

        ! Sanity check

        IF (rx > 1.0 .OR. rx < 0.0 .OR. &
            ry > 1.0 .OR. ry < 0.0 .OR. &
            rz > 1.0 .OR. rz < 0.0 ) THEN

           WRITE(*,*) '    *** Interpolation error: dvdz ***'
           WRITE(*,*) 'rx, ry, rz: ', rx, ry, rz
           WRITE(*,*) 'parcel: ', np
           WRITE(*,*) 'Aborting ...'
           CALL stopcm1

        ENDIF

       w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dvdz,pdvdz)

        ! d/dz v-tendencies

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,sgs_vz,psgs_vz)


        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,impl_vz,pimpl_vz)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,expl_vz,pexpl_vz)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,pgf_vz,ppgf_vz)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,cor_vz,pcor_vz)

!-----------------------------------------------------------------------------
! 7. dwdx
!-----------------------------------------------------------------------------

        istart = ib+1
        iend   = ie

        jstart = jb
        jend   = je

        kstart = kb
        kend   = ke+1

        DO ii = istart, iend-1
          IF ( (pdata(np,1) >= x_dwdx(ii)) .AND. (pdata(np,1) < x_dwdx(ii+1))) THEN
            i = ii
            EXIT
          ENDIF
        ENDDO

        DO jj = jstart, jend-1
          IF ( (pdata(np,2) >= y_dwdx(jj)) .AND. (pdata(np,2) < y_dwdx(jj+1)) ) THEN
            j = jj
            EXIT
          ENDIF
        ENDDO

        DO kk = kstart, kend-1
          IF ( (pdata(np,3) >= z_dwdx(kk)) .AND. (pdata(np,3) < z_dwdx(kk+1)) ) THEN
            k = kk
            EXIT
          ENDIF
        ENDDO

        rx = ( pdata(np,1)-x_dwdx(i)) * rdx
        ry = ( pdata(np,2)-y_dwdx(j)) * rdy
        rz = ( pdata(np,3)-z_dwdx(k)) / (z_dwdx(k+1) - z_dwdx(k))

        ! Sanity check

        IF (rx > 1.0 .OR. rx < 0.0 .OR. &
            ry > 1.0 .OR. ry < 0.0 .OR. &
            rz > 1.0 .OR. rz < 0.0 ) THEN

           WRITE(*,*) '    *** Interpolation error: dwdx ***'
           WRITE(*,*) 'rx, ry, rz: ', rx, ry, rz
           WRITE(*,*) 'parcel: ', np
           WRITE(*,*) 'Aborting ...'
           CALL stopcm1

        ENDIF
       w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dwdx,pdwdx)

        ! d/dx: w-tendencies

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,sgs_wx,psgs_wx)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,impl_wx,pimpl_wx)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8, expl_wx,pexpl_wx)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,pgf_wx,ppgf_wx)

!-----------------------------------------------------------------------------
! 8. dwdy
!-----------------------------------------------------------------------------

        istart = ib
        iend   = ie

        jstart = jb+1
        jend   = je

        kstart = kb
        kend   = ke+1

        DO ii = istart, iend-1
          IF ( (pdata(np,1) >= x_dwdy(ii)) .AND. (pdata(np,1) < x_dwdy(ii+1))) THEN
            i = ii
            EXIT
          ENDIF
        ENDDO

        DO jj = jstart, jend-1
          IF ( (pdata(np,2) >= y_dwdy(jj)) .AND. (pdata(np,2) < y_dwdy(jj+1)) ) THEN
            j = jj
            EXIT
          ENDIF
        ENDDO

        DO kk = kstart, kend-1
          IF ( (pdata(np,3) >= z_dwdy(kk)) .AND. (pdata(np,3) < z_dwdy(kk+1)) ) THEN
            k = kk
            EXIT
          ENDIF
        ENDDO

        rx = ( pdata(np,1)-x_dwdy(i)) * rdx
        ry = ( pdata(np,2)-y_dwdy(j)) * rdy
        rz = ( pdata(np,3)-z_dwdy(k)) / (z_dwdy(k+1) - z_dwdy(k))

       ! Sanity check

        IF (rx > 1.0 .OR. rx < 0.0 .OR. &
            ry > 1.0 .OR. ry < 0.0 .OR. &
            rz > 1.0 .OR. rz < 0.0 ) THEN

           WRITE(*,*) '    *** Interpolation error: dwdy ***'
           WRITE(*,*) 'rx, ry, rz: ', rx, ry, rz
           WRITE(*,*) 'parcel: ', np
           WRITE(*,*) 'Aborting ...'
           CALL stopcm1
        ENDIF
        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dwdy,pdwdy)

        ! d/dy: w-tendencies

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,sgs_wy,psgs_wy)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,impl_wy,pimpl_wy)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,expl_wy,pexpl_wy)

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,pgf_wy,ppgf_wy)

!-----------------------------------------------------------------------------
! 9. dwdz
!-----------------------------------------------------------------------------

        istart = ib
        iend   = ie

        jstart = jb
        jend   = je

        kstart = kb
        kend   = ke

        DO ii = istart, iend-1
          IF ( (pdata(np,1) >= x_dwdz(ii)) .AND. (pdata(np,1) < x_dwdz(ii+1)) ) THEN
            i = ii
            EXIT
          ENDIF
        ENDDO

        DO jj = jstart, jend-1
          IF ( (pdata(np,2) >= y_dwdz(jj)) .AND. (pdata(np,2) < y_dwdz(jj+1)) ) THEN
            j = jj
            EXIT
          ENDIF
        ENDDO

        DO kk = kstart, kend-1
          IF ( (pdata(np,3) >= z_dwdz(kk)) .AND. (pdata(np,3) < z_dwdz(kk+1)) ) THEN
            k = kk
            EXIT
          ENDIF
        ENDDO

        rx = ( pdata(np,1)-x_dwdz(i)) * rdx
        ry = ( pdata(np,2)-y_dwdz(j)) * rdy
        rz = ( pdata(np,3)-z_dwdz(k)) / (z_dwdz(k+1) - z_dwdz(k))

    ! Sanity check

        IF (rx > 1.0 .OR. rx < 0.0 .OR. &
            ry > 1.0 .OR. ry < 0.0 .OR. &
            rz > 1.0 .OR. rz < 0.0 ) THEN

           WRITE(*,*) '    *** Interpolation error: dwdz ***'
           WRITE(*,*) 'rx, ry, rz: ', rx, ry, rz
           WRITE(*,*) 'parcel: ', np
           WRITE(*,*) 'Aborting ...'
           CALL stopcm1

        ENDIF
        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dwdz,pdwdz)

!-----------------------------------------------------------------------------
! 10. xi
!-----------------------------------------------------------------------------

        istart = ib
        iend   = ie

        jstart = jb+1
        jend   = je

        kstart = kb+1
        kend   = ke

        DO ii = istart, iend-1
          IF ( (pdata(np,1) >= x_xi(ii)) .AND. (pdata(np,1) < x_xi(ii+1))) THEN
            i = ii
            EXIT
          ENDIF
        ENDDO

        DO jj = jstart, jend-1
          IF ( (pdata(np,2) >= y_xi(jj)) .AND. (pdata(np,2) < y_xi(jj+1)) ) THEN
            j = jj
            EXIT
          ENDIF
        ENDDO

        DO kk = kstart, kend-1
          IF ( (pdata(np,3) >= z_xi(kk)) .AND. (pdata(np,3) < z_xi(kk+1)) ) THEN
            k = kk
            EXIT
          ENDIF
        ENDDO

        rx = ( pdata(np,1)-x_xi(i)) * rdx
        ry = ( pdata(np,2)-y_xi(j)) * rdy
        rz = ( pdata(np,3)-z_xi(k)) / (z_xi(k+1) - z_xi(k))

       ! Sanity check

        IF (rx > 1.0 .OR. rx < 0.0 .OR. &
            ry > 1.0 .OR. ry < 0.0 .OR. &
            rz > 1.0 .OR. rz < 0.0 ) THEN

           WRITE(*,*) '    *** Interpolation error: xi ***'
           WRITE(*,*) 'rx, ry, rz: ', rx, ry, rz
           WRITE(*,*) 'parcel: ', np
           WRITE(*,*) 'Aborting ...'
           CALL stopcm1

        ENDIF
        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,xi,pxi)

!-----------------------------------------------------------------------------
! 11. eta
!-----------------------------------------------------------------------------

        istart = ib+1
        iend   = ie

        jstart = jb
        jend   = je

        kstart = kb+1
        kend   = ke

        DO ii = istart, iend-1
          IF ( (pdata(np,1) >= x_eta(ii)) .AND. (pdata(np,1) < x_eta(ii+1))) THEN
            i = ii
            EXIT
          ENDIF
        ENDDO

        DO jj = jstart, jend-1
          IF ( (pdata(np,2) >= y_eta(jj)) .AND. (pdata(np,2) < y_eta(jj+1)) ) THEN
            j = jj
            EXIT
          ENDIF
        ENDDO

        DO kk = kstart, kend-1
          IF ( (pdata(np,3) >= z_eta(kk)) .AND. (pdata(np,3) < z_eta(kk+1)) ) THEN
            k = kk
            EXIT
          ENDIF
        ENDDO

        rx = ( pdata(np,1)-x_eta(i)) * rdx
        ry = ( pdata(np,2)-y_eta(j)) * rdy
        rz = ( pdata(np,3)-z_eta(k)) / (z_eta(k+1) - z_eta(k))

        ! Sanity check

        IF (rx > 1.0 .OR. rx < 0.0 .OR. &
            ry > 1.0 .OR. ry < 0.0 .OR. &
            rz > 1.0 .OR. rz < 0.0 ) THEN

           WRITE(*,*) '    *** Interpolation error: eta ***'
           WRITE(*,*) 'rx, ry, rz: ', rx, ry, rz
           WRITE(*,*) 'parcel: ', np
           WRITE(*,*) 'Aborting ...'
           CALL stopcm1

        ENDIF

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,eta,peta)

!-----------------------------------------------------------------------------
! 12. zeta
!-----------------------------------------------------------------------------

        istart = ib+1
        iend   = ie

        jstart = jb+1
        jend   = je

        kstart = kb
        kend   = ke

        DO ii = istart, iend-1
          IF ( (pdata(np,1) >= x_zeta(ii)) .AND. (pdata(np,1) < x_zeta(ii+1))) THEN
            i = ii
            EXIT
          ENDIF
        ENDDO

        DO jj = jstart, jend-1
          IF ( (pdata(np,2) >= y_zeta(jj)) .AND. (pdata(np,2) < y_zeta(jj+1)) ) THEN
            j = jj
            EXIT
          ENDIF
        ENDDO

        DO kk = kstart, kend-1
          IF ( (pdata(np,3) >= z_zeta(kk)) .AND. (pdata(np,3) < z_zeta(kk+1)) ) THEN
            k = kk
            EXIT
          ENDIF
        ENDDO

        rx = ( pdata(np,1)-x_zeta(i)) * rdx
        ry = ( pdata(np,2)-y_zeta(j)) * rdy
        rz = ( pdata(np,3)-z_zeta(k)) / (z_zeta(k+1) - z_zeta(k))

        ! Sanity check

        IF (rx > 1.0 .OR. rx < 0.0 .OR. &
            ry > 1.0 .OR. ry < 0.0 .OR. &
            rz > 1.0 .OR. rz < 0.0 ) THEN

           WRITE(*,*) '    *** Interpolation error: zeta ***'
           WRITE(*,*) 'rx, ry, rz: ', rx, ry, rz
           WRITE(*,*) 'parcel: ', np
           WRITE(*,*) 'Aborting ...'
           CALL stopcm1

        ENDIF

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        CALL tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                         i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,zeta,pzeta)

!=============================================================================
! End DAHL additions
!=============================================================================

!----------------------------------------------------------------------
!  Data on u points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(np,pry).lt.yh(j) )then
          j=j-1
        endif
        if( pdata(np,prz).lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

        rx = ( pdata(np,prx)-xf(i) )/( xf(i+1)-xf(i) )
        ry = ( pdata(np,pry)-yh(j) )/( yh(j+1)-yh(j) )
        rz = ( pdata(np,prz)-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )

        ! saveit:
        rxu = rx
        rys = ry
        rzs = rz

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.0 .or. i.gt.(ni+1)   .or.        &
            j.lt.-1 .or. j.gt.(nj+1)   .or.       &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  13333b: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xf1,x3d,xf2 = ',xf(i),pdata(np,prx),xf(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),pdata(np,pry),yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),pdata(np,prz),zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni+1,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,u3d,uval)

!----------------------------------------------------------------------
!  Data on v points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(np,prx).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(np,prz).lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

        rx = ( pdata(np,prx)-xh(i) )/( xh(i+1)-xh(i) )
        ry = ( pdata(np,pry)-yf(j) )/( yf(j+1)-yf(j) )
!!!        rz = ( pdata(np,prz)-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )
        rz = rzs

        ! saveit:
        rxs = rx
        ryv = ry

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.(ni+1)   .or.       &
            j.lt.0 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  23333a: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),pdata(np,prx),xh(i+1)
          print *,'  yf1,y3d,yh2 = ',yf(j),pdata(np,pry),yf(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),pdata(np,prz),zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni,nj+1,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,v3d,vval)

!----------------------------------------------------------------------
!  Data on w points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(np,prx).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(np,pry).lt.yh(j) )then
          j=j-1
        endif

!!!        rx = ( pdata(np,prx)-xh(i) )/( xh(i+1)-xh(i) )
!!!        ry = ( pdata(np,pry)-yh(j) )/( yh(j+1)-yh(j) )
        rx = rxs
        ry = rys
        rz = ( pdata(np,prz)-zf(iflag,jflag,k) )/( zf(iflag,jflag,k+1)-zf(iflag,jflag,k) )

        ! saveit:
        rzw = rz

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.ni   .or.           &
            j.lt.-1 .or. j.gt.nj   .or.           &
            k.lt.1 .or. k.gt.nk                   )then
          print *
          print *,'  43333b: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),pdata(np,prx),xh(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),pdata(np,pry),yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zf(iflag,jflag,k),pdata(np,prz),zf(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,w3d ,wval)


! DAHL commented these sections out to avoid overwriting the parcels' locations
!      if(prkm.ge.1)then
!        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,kmh,pdata(np,prkm  ))
!        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,kmv,pdata(np,prkm+1))
!      endif
!      if(prkh.ge.1)then
!        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,khh,pdata(np,prkh  ))
!        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,khv,pdata(np,prkh+1))
!      endif
!      if( prtke.ge.1 )then
!        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,tke3d,pdata(np,prtke))
!      endif
!      if( prb.ge.1 )then
!        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,buoy,pdata(np,prb))
!      endif
!      if( prvpg.ge.1 )then
!        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,vpg(ib,jb,kb),pdata(np,prvpg))
!      endif

!----------------------------------------------------------------------
!  Data on scalar points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(np,prx).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(np,pry).lt.yh(j) )then
          j=j-1
        endif
        if( pdata(np,prz).lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

!!!        rx = ( pdata(np,prx)-xh(i) )/( xh(i+1)-xh(i) )
!!!        ry = ( pdata(np,pry)-yh(j) )/( yh(j+1)-yh(j) )
!!!        rz = ( pdata(np,prz)-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )
        rx = rxs
        ry = rys
        rz = rzs

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.ni   .or.           &
            j.lt.-1 .or. j.gt.nj   .or.           &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  15558: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *
          call stopcm1
        endif
        ENDIF

! DAHL commented out

!      if(imoist.eq.1)then
!        if(prdbz.ge.1)  &
!        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qdiag(ibdq,jbdq,kbdq,qd_dbz),pdata(np,prdbz))
!        if(prqv.ge.1)  &
!        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,q3d(ib,jb,kb,nqv),pdata(np,prqv))
!        if(prq1.ge.1)then
!          do n=nql1,nql1+(prq2-prq1)
!            call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,q3d(ib,jb,kb,n),pdata(np,prq1+(n-nql1)))
!          enddo
!        endif
!        if(prnc1.ge.1)then
!          do n=nnc1,nnc1+(prnc2-prnc1)
!            call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,q3d(ib,jb,kb,n),pdata(np,prnc1+(n-nnc1)))
!          enddo
!        endif
!        if( prqsl.ge.1 )  &
!        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dum1,pdata(np,prqsl))
!        if( prqsi.ge.1 )  &
!        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dum2,pdata(np,prqsi))
!      endif
!
!        if( prth.ge.1 )  &
!        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,th ,pdata(np,prth))
!        if( prt.ge.1 )  &
!        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,t  ,pdata(np,prt ))
!        if( prprs.ge.1 )  &
!        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,prs,pdata(np,prprs))
!   !    if( prrho.ge.1 )  &  ! DAHL
!   !     call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,rho,pdata(np,prrho))
!
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,rho,parcelrho) ! DAHL
!
!        if(prpt1.ge.1)then
!          do n=1,npt
!          call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,pt3d(ib,jb,kb,n),pdata(np,prpt1+n-1))
!          enddo
!        endif

!----------------------------------------------------------------------
!  Data on zvort points

      IF( przv.ge.1 )THEN

        i=iflag
        j=jflag
        k=kflag

        if( pdata(np,prz).lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

!!!        rx = ( pdata(np,prx)-xf(i) )/( xf(i+1)-xf(i) )
!!!        ry = ( pdata(np,pry)-yf(j) )/( yf(j+1)-yf(j) )
!!!        rz = ( pdata(np,prz)-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )
        rx = rxu
        ry = ryv
        rz = rzs

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.1 .or. i.gt.(ni+1)   .or.        &
            j.lt.1 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  15559: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *
          call stopcm1
        endif
        ENDIF

! DAHL        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,zv,pdata(np,przv))

      ENDIF

!====================================================================================
! Calculate curl of the momentum forcing functions (i.e., vorticity tendencies
! along the parcel trajectories) DAHL 20181205
!
! Method: Because of the grid stagger, we cannot not simply calculate the
! gradient field of e.g., B, V, u, w, and then interpolate.  Rather, each
! gradient is calculated on its own grid, and interpolated from its grid.  The
! interpolated values are then used to calculate the relevant expressions
! for each parcel separately.
!====================================================================================

      pxiten_sgs     = psgs_wy  - psgs_vz
      pxiten_idiff   = pimpl_wy - pimpl_vz
      pxiten_ediff   = pexpl_wy - pexpl_vz
      pxiten_pgf     = ppgf_wy  - ppgf_vz
      pxiten_cor     =          - pcor_vz  ! Tilting of f!

      petaten_sgs    = psgs_uz  - psgs_wx
      petaten_idiff  = pimpl_uz - pimpl_wx
      petaten_ediff  = pexpl_uz - pexpl_wx
      petaten_pgf    = ppgf_uz  - ppgf_wx
      petaten_cor    = pcor_uz             ! Tilting of f!

      pzetaten_sgs   = psgs_vx  - psgs_uy ! 
      pzetaten_idiff = pimpl_vx - pimpl_uy
      pzetaten_ediff = pexpl_vx - pexpl_uy
      pzetaten_pgf   = ppgf_vx  -  ppgf_uy
      pzetaten_cor   = pcor_vx  -  pcor_uy

!----------------------------------------------------------------------
!  surface variables  and  uv for parcels below lowest model level:

     ! DAHL: All assignments are now hardwired below, so don't need
     ! them here (get2d)

      IF (1 ==2) THEN
!     IF ( prznt.ge.1 .or. prust.ge.1 .or. bbc.eq.3 )THEN

        i=iflag
        j=jflag
        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( y3d.lt.yh(j) )then
          j=j-1
        endif
        call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,znt,z0)
        if( prznt.ge.1 ) pdata(np,prznt) = z0
!        if( prust.ge.1 ) call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,ust,pdata(np,prust))
      ENDIF

      IF (1 == 2) THEN
!      IF( bbc.eq.3 )THEN
        ! semi-slip lower boundary condition:
        if( z3d.lt.zh(1,1,1) )then
          ! re-calculate velocities if parcel is below lowest model level:
          !------
          ! u at lowest model level:
          i=iflag
          j=jflag
          if( y3d.lt.yh(j) )then
            j=j-1
          endif
          call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 1, 0, 1, 0,u3d(ib,jb,1),uval)
          !------
          ! v at lowest model level:
          i=iflag
          j=jflag
          if( x3d.lt.xh(i) )then
            i=i-1
          endif
          call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 1, 0, 1,v3d(ib,jb,1),vval)
          !------
          ! get u,v from (neutral) log-layer equation:
          rznt = 1.0/z0
          var = alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
          uval = uval*var
          vval = vval*var
          !------
          IF( przv.ge.1 )THEN
            do j=jflag-1,jflag+1
            do i=iflag  ,iflag+1
              z0 = 0.5*(znt(i-1,j)+znt(i,j))
              rznt = 1.0/z0
              dum3(i,j,1) = u3d(i,j,1)*alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
            enddo
            enddo
            do j=jflag  ,jflag+1
            do i=iflag-1,iflag+1
              z0 = 0.5*(znt(i,j-1)+znt(i,j))
              rznt = 1.0/z0
              dum4(i,j,1) = v3d(i,j,1)*alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
            enddo
            enddo
            do j=jflag,jflag+1
            do i=iflag,iflag+1
              dum7(i,j,1) = (dum4(i,j,1)-dum4(i-1,j,1))*rdx*uf(i)   &
                           -(dum3(i,j,1)-dum3(i,j-1,1))*rdy*vf(j)
            enddo
            enddo
            i=iflag
            j=jflag
! DAHL      call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 1, 1, 0, 0,dum7(ib,jb,1),pdata(np,przv))
          ENDIF
        endif
      ENDIF


!----------------------------------------------------------------------

! Hard-wiring all parcel assignments, and doing it at the very end of
! the routine to avoid them being overwritten by the automatic
! assignments based on namelist choices as implemented originally.

! Locations of parcels already assigned (handed over from parcel driver)
! So no need to reassign here.

pdata(np,4)  = uval
pdata(np,5)  = vval
pdata(np,6)  = wval

pdata(np,7)  = parcelrho

pdata(np,8)  = pxi
pdata(np,9)  = peta
pdata(np,10) = pzeta

pdata(np,11) = pdudx
pdata(np,12) = pdudy
pdata(np,13) = pdudz

pdata(np,14) = pdvdx
pdata(np,15) = pdvdy
pdata(np,16) = pdvdz

pdata(np,17) = pdwdx
pdata(np,18) = pdwdy
pdata(np,19) = pdwdz

pdata(np,20) = pxiten_sgs
pdata(np,21) = pxiten_idiff
pdata(np,22) = pxiten_ediff
pdata(np,23) = pxiten_pgf
pdata(np,24) = pxiten_cor

pdata(np,25) = petaten_sgs
pdata(np,26) = petaten_idiff
pdata(np,27) = petaten_ediff
pdata(np,28) = petaten_pgf
pdata(np,29) = petaten_cor

pdata(np,30) = pzetaten_sgs
pdata(np,31) = pzetaten_idiff
pdata(np,32) = pzetaten_ediff
pdata(np,33) = pzetaten_pgf
pdata(np,34) = pzetaten_cor

IF (np == -71) THEN
  print * , myid, 'PRCL VALUES', &

  pdata(np,1), &
  pdata(np,2) , &
  pdata(np,3) , &
  pdata(np,4) , &

  pdata(np,5), &
  pdata(np,6) , &
  pdata(np,7) , &

  pdata(np,8), &
  pdata(np,9), &
  pdata(np,10), &

  pdata(np,11), &
  pdata(np,12), &
  pdata(np,13), &

  pdata(np,14), &
  pdata(np,15), &
  pdata(np,16), &

  pdata(np,17), &
  pdata(np,18), &
  pdata(np,19), &
  pdata(np,20), &

  pdata(np,21), &
  pdata(np,22), &
  pdata(np,23), &
  pdata(np,24), &

  pdata(np,25), &
  pdata(np,26), &
  pdata(np,27), &
  pdata(np,28), &
  pdata(np,29), &
  pdata(np,30), &
  pdata(np,31), &
  pdata(np,32), &
  pdata(np,33), &
  pdata(np,34)

ENDIF


#ifdef MPI
      ELSE

        ! set to really small number (so we can use the allreduce command below)
        do n=4,npvals
          pdata(np,n) = -1.0e30
        enddo
#endif

      ENDIF  myprcl

    ENDDO  nploop2

!----------------------------------------------------------------------
!  communicate data
#ifdef MPI
  do np=4,npvals
    if( myid.eq.0 )then
      call MPI_REDUCE(MPI_IN_PLACE,pdata(1,np),nparcels,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    else
      call MPI_REDUCE( pdata(1,np),pdata(1,np),nparcels,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    endif
  enddo
#endif
!----------------------------------------------------------------------


      end subroutine parcel_interp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine parcel_write(prec,rtime,qname,name_prcl,desc_prcl,unit_prcl,pdata,ploc)
      use input
#ifdef NETCDF
      use writeout_nc_module, only : writepdata_nc
#endif
      implicit none

      integer, intent(inout) :: prec
      real, intent(in) :: rtime
      character(len=3), intent(in), dimension(maxq) :: qname
      character(len=40), intent(in), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl
      real, intent(in), dimension(nparcels,npvals) :: pdata
      real, intent(inout), dimension(nparcels,3) :: ploc

      integer :: n,np

!----------------------------------------------------------------------
!  write out data

    IF(myid.eq.0)THEN

      IF(output_format.eq.1)THEN
        ! GrADS format:

        string(totlen+1:totlen+22) = '_pdata.dat            '
        if(dowr) write(outfile,*) string
        open(unit=61,file=string,form='unformatted',access='direct',   &
             recl=4*npvals*nparcels,status='unknown')

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  pdata prec = ',prec

        write(61,rec=prec) ((pdata(np,n),np=1,nparcels),n=1,npvals)

!        print *, 'WRITING PARCELS:'
!DO np = 1, nparcels        
!  DO n = 1, npvals
!     print *, np, n, pdata(np,n)
!  ENDDO
!  print *, ''
!ENDDO

        close(unit=61)

#ifdef NETCDF
      ELSEIF(output_format.eq.2)THEN

        call     writepdata_nc(prec,rtime,qname,name_prcl,desc_prcl,unit_prcl,pdata,ploc(1,1))

#endif
      ENDIF
      if(dowr) write(outfile,*)

    ENDIF   ! endif for myid=0

      return
      end subroutine parcel_write


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine tri_interp(iz,jz,kz,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,s,pdata)
      use input
      implicit none

      integer :: iz,jz,kz,i,j,k
      real :: w1,w2,w3,w4,w5,w6,w7,w8
      real, dimension(1-ngxy:iz+ngxy,1-ngxy:jz+ngxy,1-ngz:kz+ngz) :: s
      real :: pdata

      pdata=s(i  ,j  ,k  )*w1    &
           +s(i+1,j  ,k  )*w2    &
           +s(i  ,j+1,k  )*w3    &
           +s(i  ,j  ,k+1)*w4    &
           +s(i+1,j  ,k+1)*w5    &
           +s(i  ,j+1,k+1)*w6    &
           +s(i+1,j+1,k  )*w7    &
           +s(i+1,j+1,k+1)*w8

      end subroutine tri_interp

!----------------------------------------------------------------------------
! Subroutine tri_interp2 // DAHL //
!----------------------------------------------------------------------------

! Changed the start/end indices to be more easily applicable to
! various fields (with different dimensions/sizes)

      subroutine tri_interp2(istart, iend, jstart, jend, kstart, kend, &
                             i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,s,pdata)
      use input

      implicit none

      integer :: istart, iend, jstart, jend, kstart, kend,i,j,k
      real :: w1,w2,w3,w4,w5,w6,w7,w8
      real, dimension(istart:iend, jstart:jend, kstart:kend) :: s
      real :: pdata

     pdata=s(i  ,j  ,k  )*w1    &
           +s(i+1,j  ,k  )*w2    &
           +s(i  ,j+1,k  )*w3    &
           +s(i  ,j  ,k+1)*w4    &
           +s(i+1,j  ,k+1)*w5    &
           +s(i  ,j+1,k+1)*w6    &
           +s(i+1,j+1,k  )*w7    &
           +s(i+1,j+1,k+1)*w8

!----------------------------------------------------------------------------
! End subroutine tri_interp2
!----------------------------------------------------------------------------

    END SUBROUTINE tri_interp2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine get2d(i,j,x3d,y3d,xh,xf,yh,yf,xs,ys,is,js,s,sval)
    use input
    implicit none

    integer, intent(in) :: i,j
    real, intent(in) :: x3d,y3d
    real, intent(in), dimension(ib:ie) :: xh
    real, intent(in), dimension(ib:ie+1) :: xf
    real, intent(in), dimension(jb:je) :: yh
    real, intent(in), dimension(jb:je+1) :: yf

    ! 0 = scalar point
    ! 1 = velocity point
    integer, intent(in) :: xs,ys
    integer, intent(in) :: is,js

    real, intent(in), dimension(ib:ie+is,jb:je+js) :: s
    real, intent(out) :: sval

    real :: wg1,wg2,wg3,wg4
    real :: x13,x23,x33,x43
    real :: w1,w2,w3,w7,rx,ry,rz

    logical, parameter :: debug = .false.

!-----------------------------------------------------------------------
      ! tri-linear interp:

      IF(xs.eq.1)THEN
        rx = ( x3d-xf(i) )/( xf(i+1)-xf(i) )
      ELSE
        rx = ( x3d-xh(i) )/( xh(i+1)-xh(i) )
      ENDIF

      IF(ys.eq.1)THEN
        ry = ( y3d-yf(j) )/( yf(j+1)-yf(j) )
      ELSE
        ry = ( y3d-yh(j) )/( yh(j+1)-yh(j) )
      ENDIF

        w1=(1.0-rx)*(1.0-ry)
        w2=rx*(1.0-ry)
        w3=(1.0-rx)*ry
        w7=rx*ry

      IF( debug )THEN
        if( rx.lt.-0.000001 .or. rx.gt.1.000001 .or.        &
            ry.lt.-0.000001 .or. ry.gt.1.000001 .or.        &
            (w1+w2+w3+w7).lt.0.999999 .or.  &
            (w1+w2+w3+w7).gt.1.000001       &
          )then
          print *,'  x3d,y3d     = ',x3d,y3d
          print *,'  i,j         = ',i,j
          print *,'  rx,ry       = ',rx,ry
          print *,'  w1,w2,w3,w7 = ',w1,w2,w3,w7
          print *,'  w1+w2+w3+w7 = ',w1+w2+w3+w7
          print *,' 22346 '
          call stopcm1
        endif
      ENDIF

      sval =s(i  ,j  )*w1    &
           +s(i+1,j  )*w2    &
           +s(i  ,j+1)*w3    &
           +s(i+1,j+1)*w7

!-----------------------------------------------------------------------

    end subroutine get2d


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine setup_parcel_vars(name_prcl,desc_prcl,unit_prcl,qname,tdef)
      use input
      implicit none

      character(len=40), intent(inout), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl
      character(len=3), intent(in), dimension(maxq) :: qname
      character(len=15), intent(inout) :: tdef

      integer :: n,n2
      character(len=8) :: text1
      character(len=30) :: text2

      prcl_out = 0

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'x'
      desc_prcl(prcl_out) = 'x position'
      unit_prcl(prcl_out) = 'm'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'y'
      desc_prcl(prcl_out) = 'y position'
      unit_prcl(prcl_out) = 'm'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'z'
      desc_prcl(prcl_out) = 'z position'
      unit_prcl(prcl_out) = 'm'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'u'
      desc_prcl(prcl_out) = 'u velocity'
      unit_prcl(prcl_out) = 'm/s'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'v'
      desc_prcl(prcl_out) = 'v velocity'
      unit_prcl(prcl_out) = 'm/s'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'w'
      desc_prcl(prcl_out) = 'w velocity'
      unit_prcl(prcl_out) = 'm/s'

      if( prth.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'th'
        desc_prcl(prcl_out) = 'potential temperature'
        unit_prcl(prcl_out) = 'K'
      endif

      if( prt.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 't'
        desc_prcl(prcl_out) = 'temperature'
        unit_prcl(prcl_out) = 'K'
      endif

      if( prprs.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'prs'
        desc_prcl(prcl_out) = 'pressure'
        unit_prcl(prcl_out) = 'Pa'
      endif

      if(prpt1.ge.1)then
        do n=1,npt
          text1='pt      '
          if(n.le.9)then
            write(text1(3:3),155) n
155         format(i1.1)
          elseif(n.le.99)then
            write(text1(3:4),154) n
154         format(i2.2)
          else
            write(text1(3:5),153) n
153         format(i3.3)
          endif

          prcl_out = prcl_out+1
          name_prcl(prcl_out) = text1
          desc_prcl(prcl_out) = 'passive tracer mixing ratio'
          unit_prcl(prcl_out) = 'kg/kg'
        enddo
      endif

      if( prqv.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'qv'
        desc_prcl(prcl_out) = 'water vapor mixing ratio'
        unit_prcl(prcl_out) = 'kg/kg'
      endif

      if(prq1.ge.1)then
        n2 = nql2
        if( iice.eq.1 ) n2 = nqs2
        do n=nql1,n2
          text1='        '
          text2='                              '
          write(text1(1:3),156) qname(n)
          write(text2(1:3),156) qname(n)
156       format(a3)

          prcl_out = prcl_out+1
          name_prcl(prcl_out) = text1
          desc_prcl(prcl_out) = text2
          unit_prcl(prcl_out) = 'kg/kg'
        enddo
      endif

      if(prnc1.ge.1)then
        do n=nnc1,nnc2
          text1='        '
          text2='                              '
          write(text1(1:3),156) qname(n)
          write(text2(1:3),156) qname(n)

          prcl_out = prcl_out+1
          name_prcl(prcl_out) = text1
          desc_prcl(prcl_out) = text2
          unit_prcl(prcl_out) = '1/kg'
        enddo
      endif

      if( prkm.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'kmh'
        desc_prcl(prcl_out) = 'horiz eddy viscosity for momentum'
        unit_prcl(prcl_out) = 'm^2/s'

        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'kmv'
        desc_prcl(prcl_out) = 'vert eddy viscosity for momentum'
        unit_prcl(prcl_out) = 'm^2/s'
      endif

      if( prkh.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'khh'
        desc_prcl(prcl_out) = 'horiz eddy diffusivity for scalars'
        unit_prcl(prcl_out) = 'm^2/s'

        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'khv'
        desc_prcl(prcl_out) = 'vert eddy diffusivity for scalars'
        unit_prcl(prcl_out) = 'm^2/s'
      endif

      if( prtke.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'tke'
        desc_prcl(prcl_out) = 'subgrid tke'
        unit_prcl(prcl_out) = 'm^2/s^2'
      endif

      if( prdbz.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'dbz'
        desc_prcl(prcl_out) = 'reflectivity'
        unit_prcl(prcl_out) = 'dBZ'
      endif

      if( prb.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'b'
        desc_prcl(prcl_out) = 'buoyancy'
        unit_prcl(prcl_out) = 'm/s/s'
      endif

      if( prvpg.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'vpg'
        desc_prcl(prcl_out) = 'vertical perturbation pressure gradient'
        unit_prcl(prcl_out) = 'm/s/s'
      endif

      if( przv.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'zvort'
        desc_prcl(prcl_out) = 'vertical vorticity'
        unit_prcl(prcl_out) = '1/s'
      endif

      if( prrho.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'rho'
        desc_prcl(prcl_out) = 'dry-air density'
        unit_prcl(prcl_out) = 'kg/m^3'
      endif

      if( prqsl.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'qsl'
        desc_prcl(prcl_out) = 'saturation mixing ratio wrt liquid'
        unit_prcl(prcl_out) = 'kg/kg'
      endif

      if( prqsi.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'qsi'
        desc_prcl(prcl_out) = 'saturation mixing ratio wrt ice'
        unit_prcl(prcl_out) = 'kg/kg'
      endif

      if( prznt.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'znt'
        desc_prcl(prcl_out) = 'surface roughness length'
        unit_prcl(prcl_out) = 'm'
      endif

      if( prust.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'ust'
        desc_prcl(prcl_out) = 'surface friction velocity'
        unit_prcl(prcl_out) = 'm/s'
      endif

!-----------------------------------------------------------------------

      if( prcl_out.gt.0 .and. output_format.eq.1 )then
        ! write GrADS descriptor file:
        call write_prclctl(name_prcl,desc_prcl,unit_prcl,tdef)
      endif

!-----------------------------------------------------------------------

      end subroutine setup_parcel_vars


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write_prclctl(name_prcl,desc_prcl,unit_prcl,tdef)
      use input
      implicit none

      character(len=40), intent(inout), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl
      character(len=15), intent(inout) :: tdef

      integer :: n,nn
      character(len=16) :: a16

      !---------------------------------------------------------------
      ! This subroutine writes the GrADS descriptor file for parcels
      !---------------------------------------------------------------

    idcheck:  &
    IF( myid.eq.0 )THEN

      string(totlen+1:totlen+22) = '_pdata.ctl            '
      if(dowr) write(outfile,*) string
      open(unit=50,file=string,status='unknown')

      sstring(baselen+1:baselen+1+12) = '_pdata.dat  '

      write(50,401) sstring
      write(50,402)
      write(50,403)
      write(50,404) nparcels
      write(50,405)
      write(50,406)
      if( prclfrq.gt.0 )then
        write(50,407) 1+int(timax/prclfrq),tdef,max(1,int(prclfrq/60.0))
      else
        write(50,407) 1000000000,tdef,max(1,int(prclfrq/60.0))
      endif

      write(50,408) prcl_out

      DO n = 1 , prcl_out
        a16 = '                '
        nn = len(trim(unit_prcl(n)))
        write(a16(2:15),214) unit_prcl(n)
        write(a16(1:1),201 )       '('
        write(a16(nn+2:nn+2),201 ) ')'
        write(50,409) name_prcl(n),desc_prcl(n),a16
      ENDDO

      write(50,410)

      close(unit=50)

    ENDIF  idcheck

201   format(a1)
214   format(a14)

401   format('dset ^',a70)
402   format('title cm1r19 output, parcel data')
403   format('undef -99999999.')
404   format('xdef ',i10,' linear 1 1')
405   format('ydef          1 linear 1 1')
406   format('zdef          1 linear 1 1')
407   format('tdef ',i10,' linear ',a15,' ',i5,'MN')
408   format('vars ',i6)
409   format(a12,' 1 99 ',a40,1x,a16)
410   format('endvars')

      end subroutine write_prclctl


!-----------------------------------------------------------------------------
! Subroutine delv
!-----------------------------------------------------------------------------

SUBROUTINE delv (ilbc, rmh, mf, xh, yh, zh, zf,                          &

                 ! Velocity gradients, velocity components, vorticity components

                 dudx, dudy, dudz, dvdx, dvdy, dvdz,                     & 
                 dwdx, dwdy, dwdz, ui, vi, w,  xi, eta, zeta,            & 

                 ! 3D accelerations
   
                 sgs_u,  sgs_v,  sgs_w,                                  &
                 impl_u, impl_v, impl_w,                                 &
                 expl_u, expl_v, expl_w,                                 &
                 pgf_u,  pgf_v,  pgf_w,                                  &
                 cor_u,  cor_v,                                          &

                 ! Antisymmetric derivatives of 3D accelerations

                 sgs_uy, impl_uy, expl_uy, pgf_uy, cor_uy,               &
                 sgs_uz, impl_uz, expl_uz, pgf_uz, cor_uz,               &
                 sgs_vx, impl_vx, expl_vx, pgf_vx, cor_vx,               &
                 sgs_vz, impl_vz, expl_vz, pgf_vz, cor_vz,               &
                 sgs_wx, impl_wx, expl_wx, pgf_wx,                       &
                 sgs_wy, impl_wy, expl_wy, pgf_wy,                       &
           
                 ! Grids associated with the various gradients

                 x_dudx, y_dudx, z_dudx, x_dudy, y_dudy, z_dudy,         &
                 x_dudz, y_dudz, z_dudz, x_dvdx, y_dvdx, z_dvdx,         &
                 x_dvdy, y_dvdy, z_dvdy, x_dvdz, y_dvdz, z_dvdz,         &
                 x_dwdx, y_dwdx, z_dwdx, x_dwdy, y_dwdy, z_dwdy,         &
                 x_dwdz, y_dwdz, z_dwdz, x_xi, y_xi, z_xi,               &
                 x_eta, y_eta, z_eta, x_zeta, y_zeta, z_zeta)

!-----------------------------------------------------------------------------
! Purpose/Method
!   Calculate the gradients of the velocity and acceleration components
!   on the staggered grid, using the centered difference across one (!)
!   grid-cell width.  This is the most accurate way to calculate centered 
!   differences, but yields values in the corners or on the faces of the grid.  
!   Hence, several different interpolation routines need to be used to ensure
!   the correct grid boxm in which the parcel resides, is identified.
!
! NOTE:
!   zh(:,:,k) is zero for k = 0.
!
! Added by J. DAHL,
! Modification history:
!
! 20190318      J. Dahl     Changed the lowest boundary conditions to one-sided
!                           to ensure that parcels within the "scalar" model
!                           domain are treated correctly.
!-----------------------------------------------------------------------------

use input
use bc_module
use comm_module

IMPLICIT NONE

! Subrourine arguments
! --------------------

INTEGER, INTENT(IN) :: &
  ilbc                         ! Flag for lower BC

REAL, INTENT(IN) ::  &
  xh(ib:ie), &
  yh(jb:je), &
  zh(ib:ie,jb:je,kb:ke), &
  zf(ib:ie,jb:je,kb:ke+1)
REAL, INTENT(IN) :: &
  ui (ib:ie+1,jb:je,   kb:ke),  &
  vi (ib:ie,  jb:je+1, kb:ke),  &
  w  (ib:ie,  jb:je,   kb:ke+1)

REAL, INTENT(IN) :: &
  sgs_u (ib:ie+1,jb:je,   kb:ke),  &
  impl_u(ib:ie+1,jb:je,   kb:ke),  &
  expl_u(ib:ie+1,jb:je,   kb:ke),  &
  pgf_u (ib:ie+1,jb:je,   kb:ke),  &
  cor_u (ib:ie+1,jb:je,   kb:ke),  &

  sgs_v (ib:ie,  jb:je+1, kb:ke),  &
  impl_v(ib:ie,  jb:je+1, kb:ke),  &
  expl_v(ib:ie,  jb:je+1, kb:ke),  &
  pgf_v (ib:ie,  jb:je+1, kb:ke),  &
  cor_v (ib:ie,  jb:je+1, kb:ke),  &

  sgs_w (ib:ie,  jb:je,   kb:ke+1), &
  impl_w(ib:ie,  jb:je,   kb:ke+1), &
  expl_w(ib:ie,  jb:je,   kb:ke+1), &
  pgf_w (ib:ie,  jb:je,   kb:ke+1)

REAL, INTENT(IN) ::      &
 rmh(ib:ie,jb:je,kb:ke), &
 mf (ib:ie,jb:je,kb:ke+1)

! Output arguments

REAL, INTENT(OUT) :: &
  dudx(ib:ie,   jb:je,   kb:ke), &
  dudy(ib:ie+1, jb+1:je, kb:ke), &
  dudz(ib:ie+1, jb:je,   kb+1:ke), &

  dvdx(ib+1:ie,jb:je+1,  kb:ke), &
  dvdy(ib:ie,  jb:je,    kb:ke), &
  dvdz(ib:ie,  jb:je+1,  kb+1:ke), &

  dwdx(ib+1:ie,jb:je,   kb:ke+1), &
  dwdy(ib:ie,  jb+1:je, kb:ke+1), &
  dwdz(ib:ie,  jb:je,   kb:ke), &

  xi  (ib:ie,   jb+1:je, kb+1:ke), &
  eta (ib+1:ie, jb:je,   kb+1:ke), &
  zeta(ib+1:ie, jb+1:je, kb:ke)

! Turbulent velocity tendency gradients

REAL, INTENT(OUT) :: &
  sgs_uy  (ib:ie+1, jb+1:je, kb:ke  ), &
  impl_uy (ib:ie+1, jb+1:je, kb:ke  ), &
  expl_uy (ib:ie+1, jb+1:je, kb:ke  ), &
  pgf_uy  (ib:ie+1, jb+1:je, kb:ke  ), &
  cor_uy  (ib:ie+1, jb+1:je, kb:ke  ), &

  sgs_uz (ib:ie+1, jb:je,   kb+1:ke), &
  impl_uz(ib:ie+1, jb:je,   kb+1:ke), &
  expl_uz(ib:ie+1, jb:je,   kb+1:ke), &
  pgf_uz (ib:ie+1, jb:je,   kb+1:ke), &
  cor_uz (ib:ie+1, jb:je,   kb+1:ke), &

  sgs_vx (ib+1:ie, jb:je+1, kb:ke  ), &
  impl_vx(ib+1:ie, jb:je+1, kb:ke  ), &
  expl_vx(ib+1:ie, jb:je+1, kb:ke  ), &
  pgf_vx (ib+1:ie, jb:je+1, kb:ke  ), &
  cor_vx (ib+1:ie, jb:je+1, kb:ke  ), &

  sgs_vz (ib:ie,   jb:je+1, kb+1:ke), & 
  impl_vz(ib:ie,   jb:je+1, kb+1:ke), &
  expl_vz(ib:ie,   jb:je+1, kb+1:ke), &
  pgf_vz (ib:ie,   jb:je+1, kb+1:ke), &
  cor_vz (ib:ie,   jb:je+1, kb+1:ke), &

  sgs_wx  (ib+1:ie, jb:je,   kb:ke+1), &
  impl_wx (ib+1:ie, jb:je,   kb:ke+1), &
  expl_wx(ib+1:ie, jb:je,   kb:ke+1), &
  pgf_wx  (ib+1:ie, jb:je,   kb:ke+1), &

  sgs_wy (ib:ie,   jb+1:je, kb:ke+1), & 
  impl_wy(ib:ie,   jb+1:je, kb:ke+1), &
  expl_wy(ib:ie,   jb+1:je, kb:ke+1), &
  pgf_wy (ib:ie,   jb+1:je, kb:ke+1)

REAL, INTENT(OUT) :: &
  x_xi(ib:ie),   &
  y_xi(jb+1:je), &
  z_xi(kb+1:ke), &

  x_eta(ib+1:ie), &
  y_eta(jb:je)  , &
  z_eta(kb+1:ke), &

  x_zeta(ib+1:ie), &
  y_zeta(jb+1:je), &
  z_zeta(kb:ke)

REAL, INTENT(OUT) :: &
  x_dudx(ib:ie), &
  y_dudx(jb:je), &
  z_dudx(kb:ke), &

  x_dudy(ib:ie+1), &
  y_dudy(jb+1:je), &
  z_dudy(kb:ke),   &

  x_dudz(ib:ie+1), &
  y_dudz(jb:je),   &
  z_dudz(kb+1:ke), &

  x_dvdx(ib+1:ie), &
  y_dvdx(jb:je+1), &
  z_dvdx (kb:ke),  &

  x_dvdy(ib:ie), &
  y_dvdy(jb:je), &
  z_dvdy(kb:ke), &

  x_dvdz(ib:ie),   &
  y_dvdz(jb:je+1), &
  z_dvdz(kb+1:ke), &

  x_dwdx(ib+1:ie), &
  y_dwdx(jb:je),   &
  z_dwdx(kb:ke+1), &

  x_dwdy(ib:ie),   &
  y_dwdy(jb+1:je), &
  z_dwdy(kb:ke+1), &

  x_dwdz(ib:ie), &
  y_dwdz(jb:je), &
  z_dwdz(kb:ke)

! Local variables
! ---------------

REAL :: &
  dudzfull(ib:ie+1, jb:je,   kb+1:ke), & ! dudz interpolated to full levels
  dvdzfull(ib:ie,   jb:je+1, kb+1:ke)    ! dvdz interpolated to full levels

REAL :: &
  dudzt(ib:ie+1,jb:je), &
  dvdzt(ib:ie,  jb:je+1)

REAL :: &
  dudz_dz2 (ib:ie+1,jb:je), &
  dutdz    (ib:ie+1,jb:je)

REAL :: &
  u(ib:ie+1,jb:je,   kb:ke),  &
  v(ib:ie,  jb:je+1, kb:ke)

REAL  ::  &
  tem

INTEGER :: &
  i, j, k


!============================================================================
! Start routine
!============================================================================

#ifdef MPI
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE)
#endif

!IF (myid == 0) print *, 'HERE'

u = ui
v = vi

! Get corner info

#ifdef MPI
!      call getcorneru3(u)
!      call bcu2(u)
!      call getcornerv3(v)
!      call bcv2(v)
!      call getcornerw3(w)
!      call bcw2(w)
#endif

! Not sure why it's applied again here.

! Apply "extrapolated" boundary condition

!IF (ilbc == 1) THEN

!  do j=jb,je
!    do i=ib,ie+1
!      tem = 0.5*(rmh(1,1,1)*mf(1,1,2))
!      u(i,j,0) = u(i,j,1)-(-1.5*u(i,j,1)+2.0*u(i,j,2)-0.5*u(i,j,3) )*tem
!    enddo
!  enddo

!  do j=jb,je+1
!    do i=ib,ie
!      tem = 0.5*(rmh(1,1,1)*mf(1,1,2))
!      v(i,j,0) =  v(i,j,1) -( -1.5*v(i,j,1) +2.0*v(i,j,2)-0.5*v(i,j,3) )*tem
!    enddo
!  enddo

!ELSEIF (ilbc == 2) THEN

! Traditional free-slip treatment (du/dz = 0)

!  do j=jb,je
!    do i=ib,ie+1
!      u(i,j,0) = u(i,j,1)
!    enddo
!  enddo

!  do j=jb,je+1
!    do i=ib,ie
!      v(i,j,0) =  v(i,j,1) 
!    enddo
!  enddo
!
!ENDIF

! Upper BCs (no parcel will ever get there, so just using the simple
! zero-gradient condition

!DO j = 0, nj+1
!  DO i = 0, ni+2
!    u(i,j,nk+1) = ui(i,j,nk)
!  ENDDO
!ENDDO

!DO j = 0, nj+2
!  DO i = 0, ni+1
!     v(i,j,nk+1) = vi(i,j,nk)
!  ENDDO
!ENDDO

!-----------------------------------------------------------------------------
! Gradients of u
!-----------------------------------------------------------------------------

! dudx
!-----

DO i = ib, ie
  DO j = jb, je
    DO k = kb, ke
      dudx(i, j, k) = (u(i+1,j,k) - u(i,j,k)) * rdx
    ENDDO
  ENDDO
ENDDO

x_dudx(ib:ie) = xh(ib:ie)
y_dudx(jb:je) = yh(jb:je)
z_dudx(kb:ke) = zh(10,10,kb:ke)

IF (myid == 43422) THEN

  DO i = ib, ie
    print '(i6, 2f12.3)', i, x_dudx(i)
  ENDDO

  print *, ''
  DO j = jb, je
    print '(i6, 2f12.3)', j,  y_dudx(j)
  ENDDO
  print *, ''

!  DO k = kb, ke
!    print '(i6, 2f12.3)', k, z_dudx(k)
!  ENDDO
ENDIF
!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!stop 888

! dudy
!-----

DO i = ib, ie+1
  DO j = jb+1, je
    DO k = kb, ke
      dudy(i,j,k) = (u(i,j,k) - u(i,j-1,k)) * rdy
    ENDDO
  ENDDO
ENDDO

! d/dy: u-tendencies

DO i = ib, ie+1
  DO j = jb+1, je
    DO k = kb, ke
      sgs_uy (i,j,k) = (sgs_u (i,j,k) - sgs_u (i,j-1,k)) * rdy
      impl_uy(i,j,k) = (impl_u(i,j,k) - impl_u(i,j-1,k)) * rdy
      expl_uy(i,j,k) = (expl_u(i,j,k) - expl_u(i,j-1,k)) * rdy
      pgf_uy (i,j,k) = (pgf_u (i,j,k) - pgf_u (i,j-1,k)) * rdy
      cor_uy (i,j,k) = (cor_u (i,j,k) - cor_u (i,j-1,k)) * rdy
    ENDDO
  ENDDO
ENDDO

x_dudy(ib)        = xh(ib) - 0.5 * dx
x_dudy(ib+1:ie+1) = xh(ib:ie) + 0.5 * dx
y_dudy(jb+1:je) = yh(jb:je-1) + 0.5 * dx
z_dudy(kb:ke)  = zh(10,10,kb:ke)

IF (myid == 24564) THEN
  DO i = ib, ie+1
    print '(i6, 2f12.3)', i, x_dudy(i)
  ENDDO

  print *, ''
  DO j = jb+1, je
    print '(i6, 2f12.3)', j,  y_dudy(j)
  ENDDO
  print *, ''

  DO k = kb, ke
    print '(i6, 2f12.3)', k, z_dudy(k)
  ENDDO
ENDIF

!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!stop 888


! du/dz
! -----

DO i = ib, ie+1
  DO j = jb, je
    DO k = kb+1, ke

      dudz(i,j,k) = (u(i,j,k) - u(i,j,k-1)) / (zh(10,10,k) - zh(10,10,k-1))

! d/dz: u-accelerations

      sgs_uz (i,j,k) = (sgs_u (i,j,k) - sgs_u (i,j,k-1)) / (zh(10,10,k) - zh(10,10,k-1))
      impl_uz(i,j,k) = (impl_u(i,j,k) - impl_u(i,j,k-1)) / (zh(10,10,k) - zh(10,10,k-1))
      expl_uz(i,j,k) = (expl_u(i,j,k) - expl_u(i,j,k-1)) / (zh(10,10,k) - zh(10,10,k-1))
      pgf_uz (i,j,k) = (pgf_u (i,j,k) - pgf_u (i,j,k-1)) / (zh(10,10,k) - zh(10,10,k-1))
      cor_uz (i,j,k) = (cor_u (i,j,k) - cor_u (i,j,k-1)) / (zh(10,10,k) - zh(10,10,k-1))

    ENDDO
  ENDDO
ENDDO

! Fix for lower boundary
! Problem: The winds (and judiag) at z = 0 are
! essentially undefined but we need dudz at dx/2, so we can safely use parcels
! that stay at/above the lowest scalar level (else, we have to discard parcels
! that descend lower than dz).  We could guess values z = 0 and calculate 
! a one-sided difference to get the dudz values at z = 0 (and then interpolate
! between z = 0 and z = dz). But, as the winds at z = 0 are ill-defined, the dudz value at
! dz/2 would be contminated by the choice of u at z = 0. Instead, we use
! the one-sided (upward) difference between dz/2 = zh(1) and zh(2), so that we get a justifiable
! dudz value at dz/2 (below that level, parcel trajectories are not robust
! anyways).  This value is the same as the centered difference at dz (dudz(i,j,kb+2)), which
! we already know.  So we'll use that value for the lower boundary:
!
! *NOTE*  
! This reasoning is based in the erroneous assumption that z_dudz(kb+1) is at
! the first half level, while it is at the bottom boundary.  The above procedure
! still makes sense but really means that we use the same vorticity value at the
! bottom boundary as at the first full level.

DO i = ib, ie+1
  DO j = jb, je
    dudz   (i,j,kb+1) = dudz   (i,j,kb+2)
    sgs_uz (i,j,kb+1) = sgs_uz (i,j,kb+2)
    impl_uz(i,j,kb+1) = impl_uz(i,j,kb+2)
    expl_uz(i,j,kb+1) = expl_uz(i,j,kb+2)
    pgf_uz (i,j,kb+1) = pgf_uz (i,j,kb+2)
    cor_uz (i,j,kb+1) = cor_uz (i,j,kb+2)
  ENDDO
ENDDO

x_dudz(ib)      = xh(ib) - 0.5 * dx
x_dudz(ib+1:ie+1) = xh(ib:ie) + 0.5 * dx
y_dudz(jb:je) = yh(jb:je)


! Using a grid that is exactlly between the full grid levels, or else the
! vertical derivatives are slightly inaccurate

z_dudz(1) = 0.0
DO k = 2, ke
  z_dudz(k) = 0.5 * (zh(10,10,k-1) + zh(10,10,k))
ENDDO

!IF (myid == 0) THEN
!  print *, 'CP1'
!  DO k = kb+1, ke
!    print *, k, zf(10,10,k),  zh(10,10,k), z_dudz(k)  
!  ENDDO
!  stop 1232
!ENDIF
 
! z_dudz(kb+1:ke)  = zf(10,10,kb+1:ke) ! Inaccurate

!  DO i = ib, ie+1
!    print '(i6, 2f12.3)', i, x_dudz(i)
!  ENDDO

!  print *, ''
!  DO j = jb, je
!    print '(i6, 2f12.3)', j,  y_dudz(j)
!  ENDDO
!  print *, ''


!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!IF (dowr) print *, 'Stopping'
!stop 888
!-----------------------------------------------------------------------------
! Gradients of v
!-----------------------------------------------------------------------------

! dvdx
! ----

DO i = ib+1, ie
  DO j = jb, je+1
    DO k = kb, ke
      dvdx(i,j,k) = (v(i,j,k) - v(i-1,j,k)) * rdx
    ENDDO
  ENDDO
ENDDO

! Turbulent tendency

DO i = ib+1, ie
  DO j = jb, je+1
    DO k = kb, ke
      sgs_vx (i,j,k) = (sgs_v (i,j,k) - sgs_v (i-1,j,k)) * rdx
      impl_vx(i,j,k) = (impl_v(i,j,k) - impl_v(i-1,j,k)) * rdx
      expl_vx(i,j,k) = (expl_v(i,j,k) - expl_v(i-1,j,k)) * rdx
      pgf_vx (i,j,k) = (pgf_v (i,j,k) - pgf_v (i-1,j,k)) * rdx
      cor_vx (i,j,k) = (cor_v (i,j,k) - cor_v (i-1,j,k)) * rdx
    ENDDO
  ENDDO
ENDDO

x_dvdx(ib+1:ie) = xh(ib+1:ie) - 0.5 * dx ! This one was wrong in earlier versions
y_dvdx(jb)        = yh(jb) - 0.5 * dy
y_dvdx(jb+1:je+1) = yh(jb:je) + 0.5 * dy
z_dvdx(kb  :ke)   = zh(10,10,kb:ke)

IF (myid == 6542) THEN
  print *, myid, 'yepyep'
  DO i = ib+1, ie
    print *, i, x_dvdx(i), xh(i)
  ENDDO
  print *, ''

  print *, ''
  DO j = jb, je+1
    print '(i6, 2f12.3)', j,  y_dvdx(j)
  ENDDO
  print *, ''

  DO k = kb, ke
    print '(i6, 2f12.3)', k, z_dvdx(k)
  ENDDO
ENDIF

!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!IF (dowr) print *, 'Stopping 112'
!stop 888

! dvdy
! ----

DO i = ib, ie
  DO j = jb, je
    DO k = kb, ke
      dvdy(i,j,k) = (v(i,j+1,k) - v(i,j,k)) * rdy
    ENDDO
  ENDDO
ENDDO

x_dvdy(ib:ie) = xh(ib:ie)
y_dvdy(jb:je) = yh(jb:je)
z_dvdy(kb:ke) = zh(10,10,kb:ke)

IF (myid == 54540) THEN
  DO i = ib, ie
    print '(i6, 2f12.3)', i, x_dvdy(i)
  ENDDO

  print *, ''
  DO j = jb, je
    print '(i6, 2f12.3)', j,  y_dvdy(j)
  ENDDO
  print *, ''

  DO k = kb, ke
    print '(i6, 2f12.3)', k, z_dvdy(k)
  ENDDO
ENDIF

!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!IF (dowr) print *, 'Stopping'
!stop 888

! dvdz
! ----

DO i = ib, ie
  DO j = jb, je+1
    DO k = kb+1, ke
      dvdz(i,j,k) = (v(i,j,k) - v(i,j,k-1)) / (zh(10,10,k) - zh(10,10,k-1))

! Turbulent and other forcings

      sgs_vz (i,j,k) = (sgs_v (i,j,k) - sgs_v (i,j,k-1)) / (zh(10,10,k) - zh(10,10,k-1))
      impl_vz(i,j,k) = (impl_v(i,j,k) - impl_v(i,j,k-1)) / (zh(10,10,k) - zh(10,10,k-1))
      expl_vz(i,j,k) = (expl_v(i,j,k) - expl_v(i,j,k-1)) / (zh(10,10,k) - zh(10,10,k-1))
      pgf_vz (i,j,k) = (pgf_v (i,j,k) - pgf_v (i,j,k-1)) / (zh(10,10,k) - zh(10,10,k-1))
      cor_vz (i,j,k) = (cor_v (i,j,k) - cor_v (i,j,k-1)) / (zh(10,10,k) - zh(10,10,k-1))

    ENDDO
  ENDDO
ENDDO

! Lower boundary conditions (see comments in the dudz section for a justification)

DO i = ib, ie
  DO j = jb, je+1
    dvdz   (i,j,kb+1) = dvdz   (i,j,kb+2)
    sgs_vz (i,j,kb+1) = sgs_vz (i,j,kb+2)
    impl_vz(i,j,kb+1) = impl_vz(i,j,kb+2)
    expl_vz(i,j,kb+1) = expl_vz(i,j,kb+2)
    pgf_vz (i,j,kb+1) = pgf_vz (i,j,kb+2)
    cor_vz (i,j,kb+1) = cor_vz (i,j,kb+2)
  ENDDO
ENDDO

x_dvdz(ib:ie) = xh(ib:ie)
y_dvdz(jb) = yh(jb) - 0.5 * dy
y_dvdz(jb+1:je+1) = yh(jb:je) + 0.5 * dy

z_dvdz(1) = 0.0
DO k = 2, ke
  z_dvdz(k) = 0.5 * (zh(10,10,k-1) + zh(10,10,k))
ENDDO

!z_dvdz(kb+1:ke)=zf(10,10,kb+1:ke) ! inaccurate (see comment in dudx section)

IF (myid == 7760) THEN
  DO i = ib, ie
    print '(i6, 2f12.3)', i, x_dvdz(i)
  ENDDO

  print *, ''
  DO j = jb, je+1
    print '(i6, 2f12.3)', j,  y_dvdz(j)
  ENDDO
  print *, ''

  DO k = kb+1, ke
    print '(i6, 2f12.3)', k, z_dvdz(k)
  ENDDO
ENDIF

!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!IF (dowr) print *, 'Stopping'
!stop 888

!-----------------------------------------------------------------------------
! Gradients of w
!-----------------------------------------------------------------------------

! dwdx
! ----

DO i = ib+1, ie
  DO j = jb, je
    DO k = kb, ke+1
      dwdx(i,j,k) = (w(i,j,k) - w(i-1,j,k)) * rdx
    ENDDO
  ENDDO
ENDDO

! Turbulent forcing

DO i = ib+1, ie
  DO j = jb, je
    DO k = kb, ke+1
      sgs_wx (i,j,k) = (sgs_w (i,j,k) - sgs_w (i-1,j,k)) * rdx
      impl_wx(i,j,k) = (impl_w(i,j,k) - impl_w(i-1,j,k)) * rdx
      expl_wx(i,j,k) = (expl_w(i,j,k) - expl_w(i-1,j,k)) * rdx
      pgf_wx (i,j,k) = (pgf_w (i,j,k) - pgf_w (i-1,j,k)) * rdx
    ENDDO
  ENDDO
ENDDO

x_dwdx(ib+1:ie) = xh(ib:ie-1) + 0.5 * dx
y_dwdx(jb:je)   = yh(jb:je)
z_dwdx(kb:ke+1) = zf(10,10,kb:ke+1)

IF (myid == 5454542) THEN
  DO i = ib+1, ie
    print '(i6, 2f12.3)', i, x_dwdx(i)
  ENDDO

  print *, ''
  DO j = jb, je
    print '(i6, 2f12.3)', j,  y_dwdx(j)
  ENDDO
  print *, ''

  DO k = kb, ke+1
    print '(i6, 2f12.3)', k, z_dwdx(k)
  ENDDO
ENDIF

!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!IF (dowr) print *, 'Stopping'
!stop 888

! dwdy
! ----

DO i = ib, ie
  DO j = jb+1, je
    DO k = kb, ke+1
      dwdy(i,j,k) = (w(i,j,k) - w(i,j-1,k)) * rdy
    ENDDO
  ENDDO
ENDDO

! Turbulent forcing

DO i = ib, ie
  DO j = jb+1, je
    DO k = kb, ke+1
      sgs_wy (i,j,k) = (sgs_w (i,j,k) - sgs_w (i,j-1,k)) * rdy
      impl_wy(i,j,k) = (impl_w(i,j,k) - impl_w(i,j-1,k)) * rdy
      expl_wy(i,j,k) = (expl_w(i,j,k) - expl_w(i,j-1,k)) * rdy
      pgf_wy (i,j,k) = (pgf_w (i,j,k) - pgf_w (i,j-1,k)) * rdy
    ENDDO
  ENDDO
ENDDO

x_dwdy(ib:ie) = xh(ib:ie)
y_dwdy(jb+1:je) = yh(jb+1:je) - 0.5 * dy
z_dwdy(kb:ke+1) = zf(10,10,kb:ke+1)

IF (myid == 8770) THEN
  DO i = ib, ie
    print '(i6, 2f12.3)', i, x_dwdy(i)
  ENDDO

  print *, ''
  DO j = jb+1, je
    print '(i6, 2f12.3)', j,  y_dwdy(j)
ENDDO
  print *, ''

  DO k = kb, ke+1
    print '(i6, 2f12.3)', k, z_dwdy(k)
  ENDDO
ENDIF

!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!IF (dowr) print *, 'Stopping'
!stop 888

! dwdz
! ----

DO i = ib, ie
  DO j = jb, je
    DO k = kb, ke
      dwdz(i,j,k) = (w(i,j,k+1) - w(i,j,k)) / (zf(i,j,k+1) - zf(i,j,k))
    ENDDO
  ENDDO
ENDDO

x_dwdz(ib:ie) = xh(ib:ie)
y_dwdz(jb:je) = yh(jb:je)
z_dwdz(kb:ke) = zh(10,10,kb:ke)

IF (myid == 76760) THEN
  DO i = ib, ie
    print '(i6, 2f12.3)', i, x_dwdz(i)
  ENDDO

  print *, ''
  DO j = jb, je
    print '(i6, 2f12.3)', j,  y_dwdz(j)
  ENDDO
  print *, ''

  DO k = kb, ke
    print '(i6, 2f12.3)', k, z_dwdz(k)
  ENDDO
ENDIF

!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!IF (dowr) print *, 'Stopping'
!stop 888

!-----------------------------------------------------------------------------
! Vorticity
!-----------------------------------------------------------------------------

! xi
! --

xi = 0.0

! Interpolating dvdz to full grid levels

DO k = 2, ke
  DO j = jb+1, je
    DO i = ib, ie
      dvdzfull(i,j,k) = dvdz(i,j,k-1) * (z_dvdz(      k)-zf    (10,10,k  ))/(z_dvdz(k)-z_dvdz(k-1)) &
                      + dvdz(i,j,k  ) * (zf    (10,10,k)-z_dvdz(      k-1))/(z_dvdz(k)-z_dvdz(k-1)) 
    ENDDO
  ENDDO
ENDDO 

dvdzfull(ib:ie,jb+1:je,kb+1) =  dvdz(ib:ie,jb+1:je,kb+1)

xi(ib:ie, jb+1:je, kb+1:ke) = dwdy(ib:ie, jb+1:je, kb+1:ke) &
                            - dvdzfull(ib:ie, jb+1:je, kb+1:ke)

x_xi(ib:ie) = xh(ib:ie)
y_xi(jb+1:je) =  yh(ib:je-1) + 0.5 * dy
z_xi (kb+1:ke)=zf(10,10,kb+1:ke)

IF (myid == 6650) THEN
  DO i = ib, ie
    print '(i6, 2f12.3)', i, x_xi(i)
  ENDDO

  print *, ''
  DO j = jb+1, je
    print '(i6, 2f12.3)', j,  y_xi(j)
  ENDDO
  print *, ''

  DO k = kb+1, ke
    print '(i6, 2f12.3)', k, z_xi(k)
  ENDDO
ENDIF

!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!IF (dowr) print *, 'Stopping t'
!stop 888

! eta
! ---

eta = 0.0

! Interpolating dvdz to full grid levels

DO k = 2, ke
  DO j = jb, je
    DO i = ib+1, ie
      dudzfull(i,j,k) = dudz(i,j,k-1) * (z_dudz(      k)-zf    (10,10,k))/(z_dudz(k)-z_dudz(k-1)) &
                      + dudz(i,j,k  ) * (zf    (10,10,k)-z_dudz(    k-1))/(z_dudz(k)-z_dudz(k-1))
    ENDDO
  ENDDO
ENDDO

dudzfull(ib+1:ie,jb:je,kb+1) =  dudz(ib+1:ie,jb:je,kb+1)

eta(ib+1:ie, jb:je, kb+1:ke) = dudzfull(ib+1:ie, jb:je, kb+1:ke) &
                             - dwdx    (ib+1:ie, jb:je, kb+1:ke)

x_eta(ib+1:ie) = xh(ib:ie-1) + 0.5 * dx
y_eta(jb:je)   = yh(jb:je)
z_eta (kb+1:ke) = zf(10,10,kb+1:ke)

IF (myid == 434326) THEN
  print *, 'BB'
  DO i = ib+1, ie
    print '(i6, 2f12.3)', i, x_eta(i)
  ENDDO

 print *, ''
  DO j = jb, je
    print '(i6, 2f12.3)', j,  y_eta(j)
  ENDDO
  print *, ''

  DO k = kb+1, ke
    print '(i6, 2f12.3)', k, z_eta(k)
  ENDDO
ENDIF

!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!IF (dowr) print *, 'Stopping '
!stop 888

! zeta
! ----

zeta = 0.0

zeta(ib+1:ie, jb+1:je, kb:ke) = dvdx(ib+1:ie, jb+1:je, kb:ke) &
                              - dudy(ib+1:ie, jb+1:je, kb:ke)

x_zeta(ib+1:ie) = xh(ib:ie-1) + 0.5 * dx
y_zeta(jb+1:je) = yh(jb:je-1) + 0.5 * dy
z_zeta(kb:ke)   = zh(10,10,kb:ke)

IF (myid == 656065) THEN
  DO i = ib+1, ie
    print '(i6, 2f12.3)', i, x_zeta(i)
  ENDDO

  print *, ''
  DO j = jb+1, je
    print '(i6, 2f12.3)', j,  y_zeta(j)
  ENDDO
  print *, ''

  DO k = kb, ke
    print '(i6, 2f12.3)', k, z_zeta(k)
  ENDDO
ENDIF

!============================================================================
! End subroutine delv
!============================================================================

END SUBROUTINE delv

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE parcel_module
