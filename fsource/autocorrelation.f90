module autocorr
    ! Module for calculating certain statistics related to spatial clustering. Depends on the
    ! spclusutils.f90 for vital statistical calcualtions, as well as other mocules used to compute
    ! various spatial statistics
    ! Ryan Martin - 2016-04-06
    ! (c) Ryan Martin 2017
    use kdtree2_module
    use rotationmatrix
    implicit none

    ! Private local variables
    logical, private :: init = .false.
    real*8, allocatable, private :: locs(:, :)
    real*8 :: rmat(3, 3)
    type(kdtree2), pointer, private :: sptree
contains
    ! initialize the autocorrelation by passing the data locations
    subroutine initautocorr( dim, ndata, dlocs, aniso )
        integer, intent(in) :: dim, ndata
        real*8, intent(in) :: dlocs(dim, ndata), aniso(5)
        !   Subroutine to initialize the kdtree for this module
        init = .false.
        ! initialize the search
        rmat = set_rmat([aniso(1), aniso(2), aniso(3)], [1.0D0, aniso(4), aniso(5)])
        ! if the tree is already constructed
        if(associated(sptree)) call kdtree2_destroy( sptree )
        ! if the locations are already allocated
        if(allocated(locs)) deallocate(locs)
        allocate(locs(dim, ndata), source=dlocs)
        ! create the KDtree
        sptree => kdtree2_create(locs, dim, .true., aniso=rmat)
        init = .true.
    end subroutine initautocorr

    ! 1. local morans  2. local getis     3. local variance
    ! 4. local mean    5. local co-variance  6. local semivariogram
    subroutine localautocorrelation( autocorrfunc, nd, nvar, rrsearch, kdtvar, z, autocorr )
        integer, intent(in)  :: nd, nvar, autocorrfunc
        logical, intent(in)  :: rrsearch   ! 1 == nn search, 2 == radius search
        real*8,  intent(in)  :: kdtvar     !variable for either nn(int) or radius(float)
        real*8,  intent(in)  :: z(nvar, nd)
        real*8,  intent(out) :: autocorr(nvar, nd)
        ! local variables
        integer :: i, j, k, nn, nf, nopt
        real*8  :: rr
        real*8, allocatable :: wts(:)
        type(kdtree2_result), allocatable, target :: dinds(:)

        if(.not. init)then
            write(*,'(a)') 'Please initialize the locations and kdtrees with initautocorr'
            return
        endif

        if(rrsearch)then
            nopt = 2
            rr = kdtvar * kdtvar
        else
            nopt = 1
            nn = int(kdtvar)
            allocate(dinds(nn),wts(nn))
        endif

        autocorr = 0.0D0
        do i = 1,nd
            if(nopt == 1)then !do the nn search
                call kdtree2_n_nearest(sptree, locs(:, i), nn, dinds)
            else
                nn = kdtree2_r_count(sptree, locs(:, i), rr)
                if(allocated(dinds)) deallocate(dinds, wts)
                if( nn < 25 )then  ! Enforce min nn for the current location:
                    nn = 25
                    allocate(dinds(nn), wts(nn))
                    call kdtree2_n_nearest(sptree, locs(:, i), nn, dinds)
                else
                    allocate(dinds(nn), wts(nn))
                    call kdtree2_r_nearest(sptree, locs(:, i), rr, nf, nn, dinds)
                endif
            endif
            ! get the weights
            wts = 0.0d0
            do j = 2, nn
                wts(j) = 1.0d0 / max( sqrt( dinds(j)%dis ), 1.0d-3 )
            enddo
            wts = wts / sum( wts )
            ! call the autocorrelation function which evaluates the similarity of the
            ! nearby points with each of the proposed metrics
            select case (autocorrfunc)
            case(1)
                autocorr(:, i) = localmorans(nn, nd, nvar, wts, z(:, i), z, dinds(:)%idx)
            case(2)
                autocorr(:, i) = localgetis(nn, nd, nvar, wts, z(:, i), z, dinds(:)%idx)
            case(3)
                autocorr(:, i) = localvariance(nn, nd, nvar, wts, z(:, i), z, dinds(:)%idx)
            case(4)
                autocorr(:, i) = localmean(nn, nd, nvar, wts, z(:, i), z, dinds(:)%idx)
            case(5)
                stop " DEPRECIATED AUTOCORRELATION FUNCTION "
                ! autocorr(:, i) = localcovariance(nn, nd, nvar, wts, z(:, i), z, dinds(:)%idx)
            case(6)
                stop " DEPRECIATED AUTOCORRELATION FUNCTION "
                ! autocorr(:, i) = localsemivariogram(nn, nd, nvar, wts, z(:, i), z, dinds(:)%idx)
            end select
        enddo

    end subroutine localautocorrelation

    function localmorans (nn, nd, nvar, wts, zi, z, idx) result ( res )
        integer, intent(in) :: nn, nd, nvar, idx(nn)
        real*8, intent(in) :: wts(nn), zi(nvar), z(nvar, nd)
        real*8  :: res(nvar), wtsum, morans_denom(nvar)
        integer :: i, j

        call weighted_variance(nn, nvar, z(:, idx), wts, morans_denom)

        ! Initialize the values
        res = 0.0D0
        wtsum = sum( wts )
        ! Loop over nearest neighbors
        do j = 1,nn
            do i = 1,nvar
                if( zi(i) < 0.0D0.and. z(i, idx(j)) < 0.0D0)then
                    res(i) = res(i) + wts(j) * ( - zi(i) * z(i, idx(j)) )
                else
                    res(i) = res(i) + wts(j) * ( zi(i) * z(i, idx(j)) )
                endif
            enddo
        enddo
        res = res / ( morans_denom * wtsum )
    end function localmorans

    ! seems like localgetis == localmean ....
    function localgetis (nn, nd, nvar, wts, zi, z, idx) result ( res )
        integer, intent(in) :: nn, nd, nvar, idx(nn)
        real*8, intent(in) :: wts(nn), zi(nvar), z(nvar, nd)
        real*8  :: res(nvar), wtsum, getis_denom(nvar)
        integer :: i, j

        ! calculate the denominator again
        getis_denom = sum( z, dim=2 )

        res = 0.0D0
        wtsum = 0.0D0
        ! Loop over nearest neighbors
        do j = 1,nn
            do i = 1,nvar
                res(i) = res(i) + wts(j) * z(i, idx(j))
                wtsum = wtsum + wts(j)
            enddo
        enddo
        res =  res / ( getis_denom * wtsum )
    end function localgetis

    function localvariance (nn, nd, nvar, wts, zi, z, idx) result ( res )
        integer, intent(in) :: nn, nd, nvar, idx(nn)
        real*8, intent(in)  :: zi(nvar), z(nvar, nd), wts(nn)
        real*8 :: res(nvar)
        ! get the weighted variance
        call weighted_variance(nn, nvar, z(:, idx), wts, res)
    end function localvariance

    function localmean (nn, nd, nvar, wts, zi, z, idx) result ( res )
        integer, intent(in) :: nn, nd, nvar, idx(nn)
        real*8, intent(in)  :: zi(nvar), z(nvar, nd), wts(nn)
        real*8 :: res(nvar)
        call weighted_mean(nn, nvar, z(:, idx), wts, res)
    end function localmean

    subroutine weighted_mean (nd, nvar, dat, wts, mean)
        integer, intent(in) :: nd, nvar
        real*8, intent(in) :: dat(nvar, nd), wts(nd)
        real*8, intent(out) :: mean(nvar)
        integer :: i
        real*8 :: wtsum
        mean = 0.0D0
        wtsum = sum(wts)
        do i = 1, nd
            mean = mean + wts(i) * dat(:, i)
        enddo
        mean = mean / wtsum
    end subroutine weighted_mean

    subroutine weighted_variance (nd, nvar, dat, wts, var)
        integer, intent(in) :: nd, nvar
        real*8, intent(in) :: dat(nvar, nd), wts(nd)
        real*8, intent(out) :: var(nvar)
        real*8 :: wtsum, mean(nvar)
        integer :: i
        ! get the weighted mean
        call weighted_mean(nd, nvar, dat, wts, mean)
        ! calc the weighted variance
        wtsum = sum( wts )
        var = 0.0D0
        do i = 1, nd
            var = var + wts(i) * ( dat(:, i) - mean ) ** 2
        enddo
        var = var / wtsum
    end subroutine weighted_variance
end module autocorr

! ! depreciated since it doesn't make much sense...
! function localcovariance (nn, nd, nvar, wts, zi, z, idx) result ( res )
!     integer, intent(in) :: nn, nd, nvar, idx(nn)
!     real*8, intent(in)  :: wts(nn), zi(nvar), z(nvar, nd)
!     real*8 :: res(nvar), z_m(nvar, nn)
!     real*8 :: var(nvar), m(nvar), mod
!     integer :: i, j, ix

!     ! get the mean and variance of the data surrounding zi
!     m = sum( z(:, idx), dim=2 ) / dble( nn )
!     var = 0.0D0
!     do i = 1, nn
!         var = var + (z(:, idx(i)) - m) ** 2
!     enddo
!     var = var / dble( nn )
!     res = 0.0D0
!     do i = 1,nn
!         do j = 1,nvar
!             mod = 1.0D0
!             if(zi(j) < 0.0D0 .and. z(j, idx(i)) < 0.0D0) mod = -1.0D0
!             res(j) = res(j) + mod * ( (zi(j)-m(j)) * (z(j, idx(i))-m(j)) ) ** 2
!         enddo
!     enddo

!     res =  var - res / dble( 2 * nd )
! end function localcovariance

! ! this one also doesnt work well from memory ...
! function localsemivariogram (nn, nd, nvar, wts, zi, z, idx) result ( res )
!     integer, intent(in) :: nn, nd, nvar, idx(nn)
!     real*8, intent(in)  :: zi(nvar), z(nvar, nd), wts(nn)
!     real*8 :: res(nvar)
!     real*8 :: m(nvar), var(nvar), wtsum
!     integer :: i, j

!     res = 0.0D0
!     do i = 1,nn
!         do j = 1,nvar
!             res(j) = res(j) + wts(i) * (zi(j) - z(j, idx(i))) ** 2
!         enddo
!     enddo
!     ! get the mean and variance of the data surrounding zi
!     call weighted_mean(nn, nvar, z(:, idx), wts, m )

!     var = 0.0D0
!     wtsum = 0.0D0
!     do i = 1, nn
!         var = var + wts(i) * (z(:, idx(i)) - m) ** 2
!         wtsum = wtsum + wts(i)
!     enddo
!     var = var / wtsum

!     res =  var - res / dble( 2 * nd )
! end function localsemivariogram

