!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine moho_stretching_honor_crust(myrank,xelm,yelm,zelm, &
                                         elem_in_crust,elem_in_mantle)

! stretching the moho according to the crust 2.0
! input:  myrank, xelm, yelm, zelm
! Dec, 30, 2009

  use constants,only: &
    NGNOD,R_EARTH_KM,R_EARTH,R_UNIT_SPHERE, &
    PI_OVER_TWO,RADIANS_TO_DEGREES,TINYVAL,SMALLVAL,ONE,USE_OLD_VERSION_5_1_5_FORMAT, &
    SUPPRESS_MOHO_STRETCHING

  use meshfem3D_par,only: &
    RMOHO_FICTITIOUS_IN_MESHER,R220,RMIDDLE_CRUST

  use meshfem3D_par,only: &
    TOPOGRAPHY

  implicit none

  integer :: myrank
  double precision,dimension(NGNOD) :: xelm,yelm,zelm
  logical,intent(inout) :: elem_in_crust,elem_in_mantle

  ! local parameters
  double precision :: r,theta,phi,lat,lon
  double precision :: vpc,vsc,rhoc,moho,elevation,gamma
  double precision :: x,y,z
  double precision :: R_moho,R_middlecrust
  integer :: ia,count_crust,count_mantle
  logical :: found_crust

  ! minimum/maximum allowed moho depths (5km/90km non-dimensionalized)
  double precision,parameter :: MOHO_MINIMUM = 5.0 / R_EARTH_KM
  ! Youyi Ruan - MARS 09/29/2015 should be able to read the maxium from outside
  !double precision,parameter :: MOHO_MAXIMUM = 90.0 / R_EARTH_KM
  double precision,parameter :: MOHO_MAXIMUM = 150.0 / R_EARTH_KM

  ! radii for stretching criteria
  R_moho = RMOHO_FICTITIOUS_IN_MESHER/R_EARTH
  R_middlecrust = RMIDDLE_CRUST/R_EARTH

  ! loops over element's anchor points
  count_crust = 0
  count_mantle = 0
  do ia = 1,NGNOD
    ! gets anchor point location
    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    if (USE_OLD_VERSION_5_1_5_FORMAT) then
      ! note: at this point, the mesh is still perfectly spherical, thus no need to
      !         convert the geocentric colatitude to a geographic colatitude
      call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
      call reduce(theta,phi)
      lat = (PI_OVER_TWO - theta) * RADIANS_TO_DEGREES
      lon = phi * RADIANS_TO_DEGREES
    else
      ! note: the moho information is given in geographic latitude/longitude (with respect to a reference ellipsoid).
      !       we need to convert the geocentric mesh positions (x/y/z) to geographic ones (lat/lon),
      !       thus correcting geocentric latitude for ellipticity
      !
      ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
      call xyz_2_rlatlon_dble(x,y,z,r,lat,lon)
    endif

    ! sets longitude bounds [-180,180]
    if (lon > 180.0d0 ) lon = lon - 360.0d0

    ! initializes
    moho = 0.d0

    ! gets smoothed moho depth
    call meshfem3D_model_crust(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,elem_in_crust)

    ! checks non-dimensionalized moho depth
    !
    ! note: flag found_crust returns .false. for points below moho,
    !          nevertheless its moho depth should be set and will be used in linear stretching
    if (moho < TINYVAL ) call exit_mpi(myrank,'Error moho depth to honor')

    if (.not. USE_OLD_VERSION_5_1_5_FORMAT) then
      ! limits moho depth to a threshold value to avoid stretching problems
      if (moho < MOHO_MINIMUM) then
        print*,'moho value exceeds minimum (in km): ',moho*R_EARTH_KM,MOHO_MINIMUM*R_EARTH_KM,'lat/lon:',lat,lon
        moho = MOHO_MINIMUM
      endif
      if (moho > MOHO_MAXIMUM) then
        print*,'moho value exceeds maximum (in km): ',moho*R_EARTH_KM,MOHO_MAXIMUM*R_EARTH_KM,'lat/lon:',lat,lon
        moho = MOHO_MAXIMUM
      endif
    endif

    ! radius of moho depth (normalized)
    moho = ONE - moho

    ! checks if moho will be honored by elements
    !
    ! note: we will honor the moho only, if the moho depth is below R_moho (~35km)
    !          or above R_middlecrust (~15km). otherwise, the moho will be "interpolated"
    !          within the element
    if (.not. SUPPRESS_MOHO_STRETCHING) then
      ! globe surface must honor topography for elements to be stretched for moho
      !
      ! note:  if no topography is honored, stretching may lead to distorted elements and invalid Jacobian
      if (TOPOGRAPHY) then
        if (moho < R_moho) then
          ! actual moho below fictitious moho
          ! elements in second layer will stretch down to honor moho topography

          elevation = moho - R_moho

          if (r >= R_moho) then
            ! point above fictitious moho
            ! gamma ranges from 0 (point at surface) to 1 (point at fictitious moho depth)
            gamma = (( R_UNIT_SPHERE - r )/( R_UNIT_SPHERE - R_moho ))
          else
            ! point below fictitious moho
            ! gamma ranges from 0 (point at R220) to 1 (point at fictitious moho depth)
            gamma = (( r - R220/R_EARTH)/( R_moho - R220/R_EARTH))

            ! since not all GLL points are exactly at R220, use a small
            ! tolerance for R220 detection, fix R220
            if (abs(gamma) < SMALLVAL) then
              gamma = 0.0d0
            endif
          endif

          if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
            call exit_MPI(myrank,'incorrect value of gamma for moho from crust 2.0')

          call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

        else  if (moho > R_middlecrust) then
          ! moho above middle crust
          ! elements in first layer will squeeze into crust above moho

          elevation = moho - R_middlecrust

          if (r > R_middlecrust) then
            ! point above middle crust
            ! gamma ranges from 0 (point at surface) to 1 (point at middle crust depth)
            gamma = (R_UNIT_SPHERE-r)/(R_UNIT_SPHERE - R_middlecrust )
          else
            ! point below middle crust
            ! gamma ranges from 0 (point at R220) to 1 (point at middle crust depth)
            gamma = (r - R220/R_EARTH)/( R_middlecrust - R220/R_EARTH )

            ! since not all GLL points are exactly at R220, use a small
            ! tolerance for R220 detection, fix R220
            if (abs(gamma) < SMALLVAL) then
              gamma = 0.0d0
            endif
          endif

          if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
            call exit_MPI(myrank,'incorrect value of gamma for moho from crust 2.0')

          call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

        endif
      endif ! TOPOGRAPHY
    endif

    ! counts corners in above moho
    ! note: uses a small tolerance
    if (r >= 0.9999d0*moho) then
      count_crust = count_crust + 1
    endif
    ! counts corners below moho
    ! again within a small tolerance
    if (r <= 1.0001d0*moho) then
      count_mantle = count_mantle + 1
    endif

  enddo

  ! sets flag when all corners are above moho
  if (count_crust == NGNOD) then
    elem_in_crust = .true.
  endif
  ! sets flag when all corners are below moho
  if (count_mantle == NGNOD) then
    elem_in_mantle = .true.
  endif

  ! small stretch check: stretching should affect only points above R220
  if (r*R_EARTH < R220) then
    print*,'Error moho stretching: ',r*R_EARTH,R220,moho*R_EARTH
    call exit_mpi(myrank,'incorrect moho stretching')
  endif

  end subroutine moho_stretching_honor_crust

!
!------------------------------------------------------------------------------------------------
!

  subroutine moho_stretching_honor_crust_reg(myrank,xelm,yelm,zelm, &
                                             elem_in_crust,elem_in_mantle)

! regional routine: for REGIONAL_MOHO_MESH adaptations
!
! uses a 3-layer crust region
!
! stretching the moho according to the crust 2.0
! input:  myrank, xelm, yelm, zelm
! Dec, 30, 2009

  use constants,only: &
    NGNOD,R_EARTH_KM,R_EARTH,R_UNIT_SPHERE, &
    PI_OVER_TWO,RADIANS_TO_DEGREES,TINYVAL,SMALLVAL,ONE,HONOR_DEEP_MOHO,USE_OLD_VERSION_5_1_5_FORMAT, &
    SUPPRESS_MOHO_STRETCHING

  use meshfem3D_par,only: &
    R220

  implicit none

  integer :: myrank

  double precision,dimension(NGNOD) :: xelm,yelm,zelm

  logical :: elem_in_crust,elem_in_mantle

  ! local parameters
  integer:: ia,count_crust,count_mantle
  double precision:: r,theta,phi,lat,lon
  double precision:: vpc,vsc,rhoc,moho
  logical:: found_crust
  double precision :: x,y,z

  ! loops over element's anchor points
  count_crust = 0
  count_mantle = 0
  do ia = 1,NGNOD

    ! anchor point location
    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    if (USE_OLD_VERSION_5_1_5_FORMAT) then
      ! note: at this point, the mesh is still perfectly spherical, thus no need to
      !         convert the geocentric colatitude to a geographic colatitude
      call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
      call reduce(theta,phi)
      lat = (PI_OVER_TWO - theta) * RADIANS_TO_DEGREES
      lon = phi * RADIANS_TO_DEGREES
    else
      ! note: the moho information is given in geographic latitude/longitude (with respect to a reference ellipsoid).
      !       we need to convert the geocentric mesh positions (x/y/z) to geographic ones (lat/lon),
      !       thus correcting geocentric latitude for ellipticity
      !
      ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
      call xyz_2_rlatlon_dble(x,y,z,r,lat,lon)
    endif

    ! sets longitude bounds [-180,180]
    if (lon > 180.d0 ) lon = lon - 360.0d0

    ! initializes
    moho = 0.d0

    ! gets smoothed moho depth
    call meshfem3D_model_crust(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,elem_in_crust)

    !=== Youyi Ruan 10/21/2016 Mars spot check
    !if (lon > 140.0 .and. lon < 143. .and. lat < -42. .and. lat > -45. .and. ia == 27) then
    !  print '(a, I4, 4F12.4)', 'crust map:', ia, lat, lon, (ONE-r)*R_EARTH_KM, moho*R_EARTH_KM 
    !endif

    if (r > 1.02d0) then
      print '(a,4F12.4)', 'DEBUG moho_stretching_honor_crust_reg(): lat, lon, r, moho =', lat, lon, r, moho
      stop 'stop for detecting large r'
    endif
    !=== end 

    ! checks moho depth
    if (abs(moho) < TINYVAL ) call exit_mpi(myrank,'Error moho depth to honor')

    moho = ONE - moho

    ! checks if moho will be honored by elements
    !
    ! note: we will honor the moho, if the moho depth is
    !         - above 15km
    !         - between 25km and 45km
    !         - below 60 km (in HONOR_DEEP_MOHO case)
    !         otherwise, the moho will be "interpolated" within the element
    if (.not. SUPPRESS_MOHO_STRETCHING) then
      ! distinguish between regions with very deep moho (e.g. Himalayan)
      if (HONOR_DEEP_MOHO) then
        call stretch_deep_moho(ia,xelm,yelm,zelm,x,y,z,r,moho)
      else
        call stretch_moho(ia,xelm,yelm,zelm,x,y,z,r,moho)
      endif
    endif

    !=== Youyi Ruan 10/21/2016 Mars spot check
    !if (lon > 140.0 .and. lon < 143. .and. lat < -42. .and. lat > -45.0 .and. ia == 27) then
    !  print '(a, I4, 4F12.4)', 'stretch:  ', ia, lat, lon, (ONE-r)*R_EARTH_KM, (ONE - moho)*R_EARTH_KM 
    !endif
    !=== end

    ! counts corners in above moho
    ! note: uses a small tolerance
    if (r >= 0.9999d0*moho) then
      count_crust = count_crust + 1
    endif
    ! counts corners below moho
    ! again within a small tolerance
    if (r <= 1.0001d0*moho) then
      count_mantle = count_mantle + 1
    endif

  enddo

  ! sets flag when all corners are above moho
  if (count_crust == NGNOD) then
    elem_in_crust = .true.
  endif
  ! sets flag when all corners are below moho
  if (count_mantle == NGNOD) then
    elem_in_mantle = .true.
  endif

  ! small stretch check: stretching should affect only points above R220
  if (r*R_EARTH < R220) then
    print*,'Error moho stretching: ',r*R_EARTH,R220,moho*R_EARTH
    call exit_mpi(myrank,'incorrect moho stretching')
  endif

  end subroutine moho_stretching_honor_crust_reg

!
!-------------------------------------------------------------------------------------------------
!

  subroutine stretch_deep_moho(ia,xelm,yelm,zelm,x,y,z,r,moho)

! honors deep moho (below 60 km), otherwise keeps the mesh boundary at r60 fixed

  use constants
  use meshfem3D_par,only: RMOHO_FICTITIOUS_IN_MESHER,R220,RMIDDLE_CRUST

  implicit none

  integer ia

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  double precision :: x,y,z

  double precision :: r,moho

  ! local parameters
  double precision :: elevation,gamma
  ! radii for stretching criteria
  !=== Youyi Ruan 10/21/2016 MARS
  double precision,parameter ::  R15=3375000.d0/R_EARTH ! 15 km  bot of 1st lyr
  double precision,parameter ::  R25=3365000.d0/R_EARTH ! 25 km   
  double precision,parameter ::  R30=3355000.d0/R_EARTH ! 35 km  bot of 2nd lyr
  double precision,parameter ::  R35=3345000.d0/R_EARTH ! 45 km
  double precision,parameter ::  R40=3335000.d0/R_EARTH ! 55 km  bot of 3rd lyr
  double precision,parameter ::  R45=3325000.d0/R_EARTH ! 65 km
  double precision,parameter ::  R50=3310000.d0/R_EARTH ! 80 km  bot of 4th lyr
  double precision,parameter ::  R55=3295000.d0/R_EARTH ! 95 km 
  double precision,parameter ::  R60=3270000.d0/R_EARTH ! 110 km moho + 10 km adjustment 

  ! === Youyi Ruan 10/21/2016 MARS
  ! checks moho position: supposed to be at 60 km

  if (RMOHO_STRETCH_ADJUSTMENT /= -10000.d0 ) &
    stop 'wrong moho stretch adjustment for stretch_deep_moho'

  if (RMOHO_FICTITIOUS_IN_MESHER/R_EARTH /= R60 ) &
    stop 'wrong moho depth '

  !=== Youyi Ruan 11/16/2016 doesn't need to check this for Mars
  ! checks middle crust position: supposed to be bottom of first layer at 15 km
  !if (RMIDDLE_CRUST/R_EARTH /= R15 ) &
  !  stop 'wrong middle crust depth'

  ! stretches mesh by moving point coordinates
  if (moho < R25 .and. moho >= R35) then
    ! moho between r25 and r35

    ! stretches mesh at R30 to moho depth
    elevation = moho - R30
    if (r >= R30 .and. r < R15) then
      gamma=((R15 - r)/(R15 - R30)) 
    else if (r < R30 .and. r > R60) then
      gamma = ((r - R60)/(R30 - R60)) ! keeps R60 fixed
      if (abs(gamma)<SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif

    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho between R25 and R35'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

  else if (moho < R35 .and. moho >= R45) then
    ! moho between R35 and r45

    ! stretch mesh at R40 to moho depth
    elevation = moho - R40
    if (r >= R40 .and. r < R15) then
      gamma=((R15 - r)/(R15 - R40))
    else if (r < R40 .and. r > R60) then
      gamma=((r - R60)/(R40 - R60)) ! keeps R60 fixed
      if (abs(gamma)<SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif

    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho between R35 and R45'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

  else if (moho < R45 .and. moho >= R55) then
    ! moho between R45 and R55

    ! stretch mesh at R50 to moho depth
    elevation = moho - R50
    if (r >= R50 .and. r < R15) then
      gamma=((R15 - r)/(R15 - R50)) 
    else if (r < R50 .and. r > R60) then
      gamma=((r - R60)/(R50 - R60)) ! keeps R60 fixed
      if (abs(gamma)<SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif

    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho between R45 and R55'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)
    
  else if (moho < R55) then
    ! moho below R55

    ! move R50 down to R55
    elevation = R55 - R50
    if (r >= R50 .and. r < R15) then 
      gamma=((R15-r)/(R15-R50))
    else if (r < R50 .and. r > R60) then
      gamma=((r - R60)/(R50 - R60)) ! keeps R60 fixed
      if (abs(gamma)<SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif

    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho from crust 2.0'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

    ! add deep moho here
    if (moho < R60) then
      ! moho below R60

      ! stretches mesh at R60 to moho
      elevation = moho - R60
      if (r < R55 .and. r >= R60) then
        gamma=(R55 - r)/(R55 - R60)
      else if (r < R60 .and. r > R220/R_EARTH) then
        gamma=(r - R220/R_EARTH)/(R60 - R220/R_EARTH)
        if (abs(gamma)<SMALLVAL) then
          gamma=0.0d0
        endif
      else
        gamma=0.0d0
      endif

      !=== Youyi Ruan MARS 10/21/2016 
      if (gamma < -0.0001d0 .or. gamma > 1.0001d0) then
        print '(a, 2F10.4,F9.6)', 'DEBUG (moho<R60)', moho*R_EARTH, r*R_EARTH, gamma
        stop 'DEBUG: moho < R60, incorrect value of gamma for moho from crust 2.0'
      endif
      !=== end 

      call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)
    endif

  else if (moho > R25) then
    ! moho above R25

    ! moves mesh at R30 up to R25
    elevation = R25 - R30
    if (r >= R30 .and. r < R15) then
      gamma=((R15 - r)/(R15 - R30)) !
    else if (r < R30 .and. r > R60) then
      gamma=(r - R60)/(R30 - R60) ! keeps r60 fixed
      if (abs(gamma)<SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif

    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho from crust 2.0'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

    ! add shallow moho here
    if (moho > R15) then
      ! moho above R15

      ! stretches mesh at R15 to moho depth
      elevation = moho - R15
      if (r >= R15) then
        gamma=(R_UNIT_SPHERE - r)/(R_UNIT_SPHERE - R15)
      else if (r < R15 .and. R > R25) then
        gamma=(r - R25)/(R15 - R25) ! keeps mesh at R25 fixed
        if (abs(gamma)<SMALLVAL) then
          gamma=0.0d0
        endif
      else
        gamma=0.0d0
      endif

      !=== Youyi Ruan MARS 10/21/2016 
      if (gamma < -0.0001d0 .or. gamma > 1.0001d0) then
        !print '(a, 3F12.4)', 'DEBUG (moho > R15)', moho*R_EARTH/1000., r*R_EARTH/1000., gamma
        !if (gamma < -0.2d0) gamma = -0.2d0
        !if (gamma > 1.02d0) gamma = 1.2d0
        stop 'DEBUG: moho > R15, incorrect value of gamma for moho above R15'
      endif
      !=== end 
      call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)
    endif
  endif

  end subroutine stretch_deep_moho

!
!-------------------------------------------------------------------------------------------------
!

  subroutine stretch_moho(ia,xelm,yelm,zelm,x,y,z,r,moho)

! honors shallow and middle depth moho, deep moho will be interpolated within elements
! mesh will get stretched down to r220

  use constants
  use meshfem3D_par,only: RMOHO_FICTITIOUS_IN_MESHER,R220,RMIDDLE_CRUST

  implicit none

  integer ia

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  double precision :: r,moho
  double precision :: x,y,z

  ! local parameters
  double precision :: elevation,gamma
  ! radii for stretching criteria
  double precision,parameter ::  R15=6356000.d0/R_EARTH
  double precision,parameter ::  R25=6346000.d0/R_EARTH
  double precision,parameter ::  R30=6341000.d0/R_EARTH
  double precision,parameter ::  R35=6336000.d0/R_EARTH
  double precision,parameter ::  R40=6331000.d0/R_EARTH
  double precision,parameter ::  R45=6326000.d0/R_EARTH
  double precision,parameter ::  R50=6321000.d0/R_EARTH
  double precision,parameter ::  R55=6316000.d0/R_EARTH
  double precision,parameter ::  R60=6311000.d0/R_EARTH


  ! checks moho position: supposed to be at 55 km
  if (RMOHO_STRETCH_ADJUSTMENT /= -15000.d0 ) &
    stop 'wrong moho stretch adjustment for stretch_moho'
  if (RMOHO_FICTITIOUS_IN_MESHER/R_EARTH /= R55 ) &
    stop 'wrong moho depth '
  ! checks middle crust position: supposed to be bottom of first layer at 15 km
  if (RMIDDLE_CRUST/R_EARTH /= R15 ) &
    stop 'wrong middle crust depth'

  ! moho between 25km and 45 km
  if (moho < R25 .and. moho > R45) then

    elevation = moho - R35
    if (r >=R35.and.r<R15) then
      gamma=((R15-r)/(R15-R35))
    else if (r<R35.and.r>R220/R_EARTH) then
      gamma = ((r-R220/R_EARTH)/(R35-R220/R_EARTH))
      if (abs(gamma)<SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif
    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho from crust 2.0'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

  else if (moho < R45) then
    ! moho below 45 km

    ! moves mesh at r35 down to r45
    elevation = R45 - R35
    if (r>= R35.and.r<R15) then
      gamma=((R15-r)/(R15-R35))
    else if (r<R35.and.r>R220/R_EARTH) then
      gamma=((r-R220/R_EARTH)/(R35-R220/R_EARTH))
      if (abs(gamma)<SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif
    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho from crust 2.0'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

  else if (moho > R25) then
    ! moho above 25km

    ! moves mesh at r35 up to r25
    elevation = R25-R35
    if (r>=R35.and.r<R15) then
      gamma=((R15-r)/(R15-R35))
    else if (r<R35.and.r>R220/R_EARTH) then
      gamma=(r-R220/R_EARTH)/(R35-R220/R_EARTH)
      if (abs(gamma)<SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif
    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho from crust 2.0'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

    ! add shallow moho here
    if (moho >R15) then
      elevation = moho-R15
      if (r>=R15) then
        gamma=(R_UNIT_SPHERE-r)/(R_UNIT_SPHERE-R15)
      else if (r<R15.and.R>R25) then
        gamma=(r-R25)/(R15-R25)
        if (abs(gamma)<SMALLVAL) then
          gamma=0.0d0
        endif
      else
        gamma=0.0d0
      endif

      call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)
    endif
  endif

  end subroutine stretch_moho

!
!-------------------------------------------------------------------------------------------------
!

  subroutine move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

! moves a point to a new location defined by gamma,elevation and r

  use constants,only: NGNOD,ONE

  implicit none

  integer ia

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  double precision :: x,y,z

  double precision :: r,elevation,gamma

  ! local parameters
  double precision :: stretch_factor

  !  stretch factor
  ! offset will be gamma * elevation
  ! scaling Cartesian coordinates xyz rather than spherical r/theta/phi involves division of offset by r
  stretch_factor = ONE + gamma * elevation/r

  ! new point location
  x = x * stretch_factor
  y = y * stretch_factor
  z = z * stretch_factor

  ! stores new point location
  xelm(ia) = x
  yelm(ia) = y
  zelm(ia) = z

  ! new radius
  r = dsqrt(xelm(ia)*xelm(ia) + yelm(ia)*yelm(ia) + zelm(ia)*zelm(ia))

  end subroutine move_point
