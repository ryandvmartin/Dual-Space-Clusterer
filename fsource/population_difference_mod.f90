! (c) Ryan Martin 2017 under MIT license
module population_difference_mod
    ! reuquires a lapack library for linking
    implicit none
    real*8, parameter, private :: PI = 3.141592653589793238462643383279502884197
contains
    ! ==============================================================================================
    ! Euclidean distance(s)
    ! ==============================================================================================
    ! Calculate the ndimensional distance given the input points pt1 and pt2
    real*8 function nd_dis(nvar, pt1, pt2)
        integer, intent(in) :: nvar
        real*8, intent(in)  :: pt1(:), pt2(:)
        nd_dis = sqrt( nd_dis_sqr( nvar, pt1, pt2 ) )
    end function nd_dis

    ! Calculate the ndimensional squared distance given the input points pt1 and pt2
    real*8 function nd_dis_sqr(nvar, pt1, pt2)
        integer, intent(in) :: nvar
        real*8, intent(in)  :: pt1(:), pt2(:)
        real*8 :: tsum
        integer :: i
        tsum = 0.0d0
        do i = 1,nvar
            tsum = tsum + ( pt2(i) - pt1(i) ) ** 2
        enddo
        nd_dis_sqr = tsum
    end function nd_dis_sqr

    ! get the KLD by computing the probabilities for:
    !   x under x, x under y, y under y and y under x
    function kl_distance( nvar, n, m, x, y, band ) result ( distance )
        integer, intent(in) :: nvar, n, m
        real*8, intent(in) :: x(nvar, n), y(nvar, m), band
        real*8 :: xt(nvar, n), yt(nvar, m), u_x(nvar), u_y(nvar), &
                  cov_x(nvar, nvar), cov_y(nvar, nvar), inv_x(nvar, nvar), inv_y(nvar, nvar)
        integer :: ii, fstat, i, errcode1, errcode2
        real*8 :: distance, px(n), qx(n), py(m), qy(m)
        real*8, parameter :: regconst=1.0d-4

        ! get the mean and the centered variables
        u_x = sum( x, dim=2 ) / dble( n )
        u_y = sum( y, dim=2 ) / dble( m )
        do ii = 1, n
            xt(:, ii) = x(:, ii) - u_x
        enddo
        do ii = 1, m
            yt(:, ii) = y(:, ii) - u_y
        enddo
        ! calculate the covariance matrices
        cov_x = matmul( xt, transpose(xt) ) / dble( n )
        cov_y = matmul( yt, transpose(yt) ) / dble( m )
        ! add the regularization constant to the diagonal of the covariance matrix
        do i = 1, nvar
            cov_x(i, i) = cov_x(i, i) + regconst
            cov_y(i, i) = cov_y(i, i) + regconst
        enddo
        inv_x = cov_x
        inv_y = cov_y

        ! invert the covariance matrices
        call invert( inv_x, nvar, errcode1 )
        call invert( inv_y, nvar, errcode2 )

        ! scale the covariance by the bandwidth
        inv_x = inv_x / ( 2.0d0 * band )
        inv_y = inv_y / ( 2.0d0 * band )
        do i = 1, n
            ! density for population x, from x
            px(i) = gauss_density_at_point( nvar, n, x(:, i), x, cov_x, inv_x )
            ! density for population x, from y
            qx(i) = gauss_density_at_point( nvar, n, x(:, i), y, cov_y, inv_y )
        enddo
        do i = 1, m
            ! density for population y, from y
            py(i) = gauss_density_at_point( nvar, m, y(:, i), y, cov_y, inv_y )
            ! density for population y, from x
            qy(i) = gauss_density_at_point( nvar, m, y(:, i), x, cov_x, inv_x )
        enddo

        ! take the average of the 'forward' and 'reverse' KLD
        distance = 0.5d0 + ( sum( px * log( px / qx ) ) + sum( py * log( py / qy ) ) )

    end function kl_distance

    function gauss_density_at_point( nvar, n, pt, pop, cov, invcov ) result ( dens )
        integer, intent(in) :: nvar, n
        real*8, intent(in) :: pt(nvar), pop(nvar, n), cov(nvar, nvar), invcov(nvar, nvar)
        real*8 :: dens, dist(nvar), norm
        integer :: i
        dens = 0.0d0
        do i = 1, n
            dist = pop(:, i) - pt
            dens = dens + exp( - mdist_sqr(nvar, dist, invcov) )
        enddo
        norm = 0.0D0
        do i = 1, nvar
            norm = norm + sum( cov(:, i) ** 2 )
        enddo
        norm = sqrt(norm)
        norm = 1.0d0 / ( (2.0d0*PI)**(dble(nvar)/2.0d0) * norm )
        dens = max(1d-12, norm * dens)
    end function gauss_density_at_point

    ! Wards distance, the increase in variance obtained by combining populations x and y
    function wards_distance( nvar, n, m, x, y ) result ( distance )
        integer, intent(in) :: nvar, n, m
        real*8, intent(in) :: x(nvar, n), y(nvar, m)
        real*8 :: SSE_xy, SSE_x, SSE_y, xm(nvar), ym(nvar), xym(nvar), distance
        integer :: i
        ! get the mean of each population
        xm = sum( x, dim=2 ) / dble( n )
        ym = sum( y, dim=2 ) / dble( m )
        xym = 0.5d0 * ( xm + ym )

        ! initialize the SSE variables
        SSE_xy = 0.0D0
        SSE_x = 0.0D0
        SSE_y = 0.0D0
        ! calculate the SSE's
        do i = 1, n
            SSE_xy = SSE_xy + sum( (x(:, i) - xym)**2 )
            SSE_x  = SSE_x  + sum( (x(:, i) - xm )**2 )
        enddo
        do i = 1, m
            SSE_xy = SSE_xy + sum( (y(:, i) - xym)**2 )
            SSE_y  = SSE_y  + sum( (y(:, i) - ym )**2 )
        enddo

        ! finalize the statistic
        distance = (SSE_xy - ( SSE_x + SSE_y )) / SSE_xy
    end function wards_distance

    ! energy distance from Rizzo 2016
    ! assume same number of variables (for now), and n and m are the size of x and y respectively
    function energy_statistic( nvar, n, m, x, y ) result ( distance )
        integer, intent(in) :: nvar, n, m
        real*8, intent(in) :: x(nvar, n), y(nvar, m)
        real*8 :: a, b, c, xt(nvar), yt(nvar), distance
        integer :: i, j

        a = 0.0D0
        b = 0.0D0
        c = 0.0D0
        ! calculate a
        do i = 1, n
            xt = x(:, i)
            do j = 1, m
                a = a + nd_dis( nvar, xt, y(:, j) )
            enddo
        enddo
        a = a / dble( n * m )

        ! calculate b
        do i = 1, n
            xt = x(:, i)
            do j = 1, n
                b = b + nd_dis( nvar, xt, x(:, j) )
            enddo
        enddo
        b = b / dble( n**2 )

        ! calculate c
        do i = 1, m
            yt = y(:, i)
            do j = 1, m
                c = c + nd_dis( nvar, yt, y(:, j) )
            enddo
        enddo
        c = c / dble( m**2 )

        ! the final distance metric
        distance = (2.0D0 * a - b - c) / ( 2.0D0 * a )
    end function energy_statistic

    ! Modified energy distance from Rizzo 2016 for the mahalanobis distance
    ! assume same number of variables (for now), and n and m are the size of x and y respectively
    function energy_statistic_mdist( nvar, n, m, x, y, regconst ) result ( distance )
        integer, intent(in) :: nvar, n, m
        real*8, intent(in) :: x(nvar, n), y(nvar, m), regconst
        real*8 :: xt(nvar, n), yt(nvar, m), xty(nvar, n), ytx(nvar, m), &
                  u_x(nvar), u_y(nvar), cov_x(nvar, nvar), cov_y(nvar, nvar)
        real*8 :: a, b, c, distance
        integer :: ii, fstat, i, errcode1, errcode2

        ! get the mean and generate the temporary x-x, x-y, y-y, y-x vectors
        u_x = sum( x, dim=2 ) / size( x, 2 )
        u_y = sum( y, dim=2 ) / size( y, 2 )
        do ii = 1, n
            xt(:, ii) = x(:, ii) - u_x
            xty(:, ii) = x(:, ii) - u_y
        enddo
        do ii = 1, m
            yt(:, ii) = y(:, ii) - u_y
            ytx(:, ii) = y(:, ii) - u_x
        enddo
        cov_x = matmul( xt, transpose(xt) ) / dble( n )
        cov_y = matmul( yt, transpose(yt) ) / dble( m )

        ! add the regularization constant to the diagonal of the covariance matrix
        do i = 1, nvar
            cov_x(i, i) = cov_x(i, i) + regconst
            cov_y(i, i) = cov_y(i, i) + regconst
        enddo

        ! invert the covariance matrices
        call invert( cov_x, nvar, errcode1 )
        call invert( cov_y, nvar, errcode2 )

        ! if(errcode1 /= 0 .or. errcode2 /= 0) stop " ERROR inverting covariance matrix "

        ! calculate A by getting the mahalanobis distance from x to y, then y to x
        a = 0.0D0
        ! get (y-u_x)T Cx-1 (y-u_x)
        do ii = 1, m
            a = a + mdist( nvar, ytx(:, ii), cov_x )
        enddo
        ! get (x-u_y)T Cy-1 (x-u_y)
        do ii = 1, n
            a = a + mdist( nvar, xty(:, ii), cov_y )
        enddo
        a = a / dble( n * m )

        ! get the distance sum betweenx and x
        b = 0.0D0
        do ii = 1, n
            b = b + mdist( nvar, xt(:, ii), cov_x )
        enddo
        b = b / dble( n**2 )

        ! get the distance between y and y
        c = 0.0D0
        do ii = 1, m
            c = c + mdist( nvar, yt(:, ii), cov_y )
        enddo
        c = c / dble( m**2 )

        distance = (a - b - c) / ( a )
    end function energy_statistic_mdist

    !Mahalanobis distance using the inverted covariance matrix
    function mdist( d, x, cinv ) result ( mdistance )
        integer, intent(in) :: d
        real*8, intent(in) :: x(d, 1)    ! The centered value i.e. = x - u_x
        real*8, intent(in) :: cinv(d, d) ! the inverted covariance matrix
        real*8 :: x_t(1, d), distance(1, 1), mdistance
        x_t(1, :) = x(:, 1)
        distance = matmul( x_t, matmul( cinv, x ) )
        mdistance = sqrt( abs(distance(1, 1)) )
    end function mdist

    !Mahalanobis distance using the inverted covariance matrix
    function mdist_sqr( d, x, cinv ) result ( mdistance )
        integer, intent(in) :: d
        real*8, intent(in) :: x(d, 1)    ! The centered value i.e. = x - u_x
        real*8, intent(in) :: cinv(d, d) ! the inverted covariance matrix
        real*8 :: x_t(1, d), distance(1, 1), mdistance
        x_t(1, :) = x(:, 1)
        distance = matmul( x_t, matmul( cinv, x ) )
        mdistance = abs(distance(1, 1))
    end function mdist_sqr

    ! subroutine to invert matrix `A` inplace, requires linking against a lapack library
    subroutine invert ( A, n, info )
        integer, intent (in) :: n
        real*8, intent (inout) :: A(n, n)
        integer, intent(out) :: info
        integer :: ipiv(n), LDA, LWORK
        real*8, allocatable :: WORK(:)
        LDA = n
        LWORK = 3*n-1
        allocate( WORK(LWORK), source=0.0d0)
        info = 0
        call DGETRF( n, n, A, LDA, ipiv, info )
        call DGETRI( n, A, LDA, ipiv, WORK, LWORK, info )
    end subroutine invert

end module population_difference_mod
