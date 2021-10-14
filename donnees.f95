! =====================================================================
! Module contenant les fonctions :
!
! f
! u
! phi
! subd
! base
! phiP2
! =====================================================================

MODULE donnees

    use math

    implicit none

    real(rp), dimension(2), parameter :: d_phi1 = (/-1.0_rp, -1.0_rp/)
    real(rp), dimension(2), parameter :: d_phi2 = (/ 1.0_rp,  0.0_rp/)
    real(rp), dimension(2), parameter :: d_phi3 = (/ 0.0_rp,  1.0_rp/)

contains

    ! -------------------------------------------------------------------------------------------------------
    ! -------------------------------------------------------------------------------------------------------
    function f(x, y)
        implicit none

        ! PARAMÈTRES
        real(rp), intent(in) :: x
        real(rp), intent(in) :: y

        ! VARIABLE RETOURNÉE
        real(rp) :: f

        ! dirichlet homogène
        !f = 2 * x * (10 - x) + 2 * y * (10 - y)

        ! dirichlet non homogène
        f = cos(x) + sin(y)

        ! cas affine
        !f = 0.0_rp
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! condition au bord
    ! -------------------------------------------------------------------------------------------------------
    function t0(x, y)
        implicit none

        ! paramètres
        real(rp), intent(in) :: x, y

        ! return
        real(rp) :: t0

        t0 = uExacte(x, y)
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! solution exacte pour f(x, y) = 2x(10 - x) + 2y(10 - y)
    ! -------------------------------------------------------------------------------------------------------
    function uExacte(x, y)
        implicit none

        ! paramètres
        real(rp), intent(in) :: x, y

        ! return
        real(rp) :: uExacte

        ! dirichlet homogène
        !uExacte = x * (10 - x) * y * (10 - y)

        ! dirichlet non homogène
        uExacte = cos(x) + sin(y)

        ! cas affine
        !uExacte = x + y
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! fonctions de base locales sur elt fini de référence (0, 0), (1, 0), (0, 1)
    ! -------------------------------------------------------------------------------------------------------
    function phi1(x, y)
        implicit none

        ! paramètres
        real(rp) :: x, y

        ! return
        real(rp) :: phi1

        phi1 = 1 - x - y
    end function

    function phi2(x, y)
        implicit none

        ! paramètres
        real(rp) :: x, y

        ! return
        real(rp) :: phi2

        phi2 = x
    end function

    function phi3(x, y)
        implicit none

        ! paramètres
        real(rp) :: x, y

        ! return
        real(rp) :: phi3

        phi3 = y
    end function

END MODULE donnees
