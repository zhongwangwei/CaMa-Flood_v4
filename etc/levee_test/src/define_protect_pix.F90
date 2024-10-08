      program downscale_flddph
! ===============================================
      implicit none
! CaMa-Flood parameters       
      character*256            ::  param                         !! river map parameters
      integer                  ::  iXX, iYY
      integer                  ::  nXX, nYY                      !! grid number (river network map)
      integer                  ::  nflp                          !! floodplain layers
      real*8                   ::  gsize                         !! grid size [deg]
      real*8                   ::  west, east, north, south      !! domain (river network map)

! HydroSHEDS parameters                                        !! (from location.txt)
      integer                  ::  i, narea               !! area ID
      character*256            ::  area                          !! area code
      integer                  ::  ix, iy, jx, jy                        
      integer                  ::  nx, ny                        !! grid number (hires data)
      real*8                   ::  csize                         !! size of pixel [deg]
      real*8                   ::  lon_ori                       !! west  edge
      real*8                   ::  lon_end                       !! east  edge
      real*8                   ::  lat_ori                       !! north edge
      real*8                   ::  lat_end                       !! south edge

      real*8                   ::  west2, east2, north2, south2      !! output domain 
      integer                  ::  mx, my

      character*256            ::  list_loc
!
      integer*4,allocatable    ::  nextXX(:,:)                    !! downstream (jXX,jYY)
      integer*2,allocatable    ::  catmXX(:,:), catmYY(:,:)        !! catchment (iXX,iYY) of pixel (ix,iy)
      integer*1,allocatable    ::  catmZZ(:,:)                    !! catchment Z layer (100 levels)

      real,allocatable         ::  flddif(:,:), rivwth(:,:), hand(:,:)      !! height above channel [m]
      real,allocatable         ::  lon(:), lat(:)

      real,allocatable         ::  levfrc(:,:)                   !! levee unprotected fraction
!
      real,allocatable         ::  protect(:,:)                    !! downscaled flood depth [m]
      real,allocatable         ::  flddif2(:,:),hand2(:,:)                    !! downscaled flood depth [m]

!
      character*256            ::  mapdir, hires
      parameter                   (mapdir='./map/')              !! map directory (please make a symbolic link)
      character*256            ::  fnextxy, fflood, rfile
      character*256            ::  flevfrc
      character*256            ::  buf
      integer                  ::  ios
! ===============================================
! downscale target domain
      call getarg(1,buf)
        read(buf,*) west2
      call getarg(2,buf)
        read(buf,*) east2
      call getarg(3,buf)
        read(buf,*) south2
      call getarg(4,buf)
        read(buf,*) north2

      call getarg(5,hires)       !! downscale resolution

      call getarg(6,flevfrc)       !! downscale resolution

      param=trim(mapdir)//'params.txt'

      open(11,file=param,form='formatted')
      read(11,*) nXX
      read(11,*) nYY
      read(11,*) nflp
      read(11,*) gsize
      read(11,*) west
      read(11,*) east
      read(11,*) south
      read(11,*) north
      close(11)

! CaMa-Flood simulation domain
      east =west +real(nXX)*gsize
      south=north-real(nXX)*gsize

      if( trim(hires)=='1sec' )then
        csize=1./3600.
      elseif( trim(hires)=='3sec' )then
        csize=1./1200.
      elseif( trim(hires)=='5sec' )then
        csize=1./720.
      elseif( trim(hires)=='15sec' )then
        csize=1./240.
      elseif( trim(hires)=='30sec' )then
        csize=1./120.
      elseif( trim(hires)=='1min' )then
        csize=1./60.
      else
        stop
      endif

! calculate number of cell size
      mx=nint( (east2 -west2 )/csize )   !! downscale domain x*y
      my=nint( (north2-south2)/csize )

      if( mx*my>525000000 )then
        print *, 'downscale domain too large: mx*my*4>integer limit'
        stop
      endif

      print *, 'domain:', west2, east2, south2, north2, mx, my

! ==========

      allocate(nextXX(nXX,nYY),levfrc(nXX,nYY))
      allocate(protect(mx,my),flddif2(mx,my),hand2(mx,my))
      protect(:,:)=-9999

! ===============================================
      fnextxy=trim(mapdir)//'nextxy.bin'

      open(11, file=fnextxy, form='unformatted', access='direct', recl=4*nXX*nYY)
      read(11,rec=1) nextXX
      close(11)

      open(12, file=flevfrc, form='unformatted', access='direct', recl=4*nXX*nYY)
      read(12,rec=1) levfrc
      close(12)


! open hires files
      list_loc=trim(mapdir)//trim(hires)//'/location.txt'
      open(11,file=list_loc,form='formatted')
      read(11,*) narea
      read(11,*)

      do i=1, narea
        read(11,*) buf, area, lon_ori, lon_end, lat_end, lat_ori, nx, ny, csize

        if( lon_end<west2 .or. lon_ori>east2 .or.lat_ori<south2 .or. lat_end>north2 ) goto 1090  !! out of domain
        allocate(catmXX(nx,ny),catmYY(nx,ny),catmZZ(nx,ny),flddif(nx,ny),rivwth(nx,ny),hand(nx,ny))
        allocate(lon(nx),lat(ny))
  
        rfile=trim(mapdir)//trim(hires)//'/'//trim(area)//'.catmxy.bin'
        print *, rfile
        open(21,file=rfile,form='unformatted',access='direct',recl=2*nx*ny,status='old',iostat=ios)
        if( ios==0 )then
          read(21,rec=1) catmXX
          read(21,rec=2) catmYY
          close(21)
        else
          print *, '*******************'
          print *, 'no data: ', rfile
          stop
        endif

        rfile=trim(mapdir)//trim(hires)//'/'//trim(area)//'.catmz100.bin'
        print *, rfile
        open(21,file=rfile,form='unformatted',access='direct',recl=1*nx*ny,status='old',iostat=ios)
        if( ios==0 )then
          read(21,rec=1) catmZZ
          close(21)
        else
          print *, '*******************'
          print *, 'no data: ', rfile
          stop
        endif

        rfile=trim(mapdir)//trim(hires)//'/'//trim(area)//'.flddif.bin'
        open(21,file=rfile,form='unformatted',access='direct',recl=4*nx*ny,status='old',iostat=ios)
        if( ios==0 )then
          read(21,rec=1) flddif
          close(21)
        else
          print *, '*******************'
          print *, 'no data: ', rfile
          stop
        endif

  
        rfile=trim(mapdir)//trim(hires)//'/'//trim(area)//'.rivwth.bin'
        open(21,file=rfile,form='unformatted',access='direct',recl=4*nx*ny,status='old',iostat=ios)
        if( ios==0 )then
          read(21,rec=1) rivwth
          close(21)
        else
          print *, '*******************'
          print *, 'no data: ', rfile
          stop
        endif

        rfile=trim(mapdir)//trim(hires)//'/'//trim(area)//'.hand.bin'
        open(21,file=rfile,form='unformatted',access='direct',recl=4*nx*ny,status='old',iostat=ios)
        if( ios==0 )then
          read(21,rec=1) hand
          close(21)
        else
          print *, '*******************'
          print *, 'no data: ', rfile
          stop
        endif

        do ix=1, nx
          lon(ix)=lon_ori+(real(ix)-0.5)*csize
          if( lon(ix)>=180. ) lon(ix)=lon(ix)-360.
          if( lon(ix)<-180. ) lon(ix)=lon(ix)+360.
        end do
        do iy=1, ny
          lat(iy) =lat_ori-(real(iy)-0.5)*csize
        end do

        do iy=1, ny
          do ix=1, nx
            if( lon(ix)>west2 .and. lon(ix)<east2 .and. lat(iy)>south2 .and. lat(iy)<north2 )then
              jx=int( (lon(ix)-west2 )/csize )+1
              jy=int( (north2-lat(iy))/csize )+1

              flddif2(jx,jy)=flddif(ix,iy)

              if( hand(ix,iy)>0 )then
                hand2(jx,jy)  =hand(ix,iy)**0.3
                hand2(jx,jy)  =min(hand2(jx,jy),5.)
              else
                hand2(jx,jy)=-1
              endif

              if( catmXX(ix,iy)>0 )then
                protect(jx,jy)=0                                                    !! 
                iXX=catmXX(ix,iy)
                iYY=catmYY(ix,iy)
                if( catmZZ(ix,iy)>levfrc(iXX,iYY)*100 )then  !! protected pixel
                  protect(jx,jy)=5
                endif

                if( rivwth(ix,iy)/=-9999 .and. rivwth(ix,iy)/=0 )then !! permanent water
                  if( protect(jx,jy)==0 ) protect(jx,jy)=1
                endif
              elseif( catmXX(ix,iy)/=-9999 )then
                protect(jx,jy)=-1
              endif
            endif
          end do
        end do

        deallocate(catmXX,catmYY,catmZZ,flddif,rivwth,hand)
        deallocate(lon,lat)

 1090   continue
      enddo

      hand2(1,1)= 5
      hand2(1,2)=-1

      fflood='./data/protect_pix.bin'
      open(11, file=fflood, form='unformatted', access='direct', recl=4*mx*my)
      write(11,rec=1) protect
      write(11,rec=2) flddif2
      close(11)

      fflood='./data/hand.bin'
      open(11, file=fflood, form='unformatted', access='direct', recl=4*mx*my)
      write(11,rec=1) hand2
      close(11)

      end program downscale_flddph

