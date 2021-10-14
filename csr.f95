!=====================================================================
! Module contenant les fonctions :
!
! constructCsr
! getCsr
! affCsr
! prodCsrVec
! gaussSeidelCsr
! =====================================================================

MODULE csr

    use math

    implicit none


    type :: CsrMat
        integer :: n, m, nCoeff
        real(rp), dimension(:), allocatable :: tmat
        integer, dimension(:), allocatable :: jposi, jvcell
    end type CsrMat

contains

    ! ---------------------------------------------------------------------
    ! constructeur d'une variable ce type CsrMat (matrice creuse)
    ! ---------------------------------------------------------------------
    ! n : nombre de lignes
    ! m : nombre de colonnes
    ! nCoeff : nombre de valeurs non nulles
    ! tmat, jposi, jvcell : caractérise la matrice (voir stockage csr)
    ! mat : matrice creuse retournée
    subroutine constructCsr(n, m, nCoeff, tmat, jposi, jvcell, mat)
        implicit none

        integer, intent(in) :: n, m, nCoeff
        real(rp), dimension(nCoeff), intent(in) :: tmat
        integer, dimension(n + 1), intent(in) :: jposi
        integer, dimension(nCoeff), intent(in) :: jvcell
        type(CsrMat), intent(out) :: mat

        mat = CsrMat(n, m, nCoeff, tmat, jposi, jvcell)
    end subroutine constructCsr

    subroutine rmCsr(mat)
        implicit none

        ! paramètres
        type(CsrMat), intent(inout) :: mat

        if (allocated(mat%tmat)) then
            deallocate(mat%tmat)
        end if

        if (allocated(mat%jposi)) then
            deallocate(mat%jposi)
        end if

        if (allocated(mat%jvcell)) then
            deallocate(mat%jvcell)
        end if
    end subroutine



    ! ---------------------------------------------------------------------
    ! fournit la valeur d'un elt de mat à la position (i, j)
    ! ---------------------------------------------------------------------
    ! mat : matrice creuse
    ! i, j : indices de ligne et colonne
    ! elt : valeur retournée
    subroutine getCsr(mat, i, j, elt)
        implicit none

        type(CsrMat), intent(in) :: mat
        integer, intent(in) :: i, j
        real(rp), intent(out) :: elt

        ! variables locales
        integer :: k

        elt = 0.d0

        do k = mat%jposi(i), mat%jposi(i + 1) - 1
            if (mat%jvcell(k) == j) then
                elt = mat%tmat(k)
                exit
            end if
        end do
    end subroutine getCsr



    ! ---------------------------------------------------------------------
    ! affiche une matrice creuse en console
    ! ---------------------------------------------------------------------
    ! mat : matrice creuse
    subroutine affCsr(mat)
        implicit none

        ! paramètres
        type(CsrMat), intent(in) :: mat

        ! variables locales
        real(rp), dimension(mat%m) :: temp
        integer :: i, j

        do i = 1, mat%n
            do j = 1, mat%m
                call getCsr(mat, i, j, temp(j))
            end do
            print "(100F7.2)", temp
        end do
    end subroutine affCsr



    ! ---------------------------------------------------------------------
    ! produit matrice-creuse vecteur
    ! ---------------------------------------------------------------------
    ! mat : matrice creuse de taille n * m
    ! vec : vecteur de taille m
    ! prod : vecteur retourné de taille n
    subroutine prodCsrVec(mat, vec, prod)
        implicit none

        ! paramètres
        type(CsrMat), intent(in) :: mat
        real(rp), dimension(mat%m), intent(in) :: vec
        real(rp), dimension(mat%n), intent(out) :: prod

        ! variables locales
        integer :: i, j, k
        real(rp) :: temp1, temp2

        ! coordonnée j de prod
        do i = 1, mat%n

            temp2 = 0.d0

            ! produit d'une ligne et du vecteur
            !do j = mat%jvcell(mat%jposi(i)), mat%jvcell(mat%jposi(i + 1) - 1)
            do k = mat%jposi(i), mat%jposi(i + 1) - 1
                j = mat%jvcell(k)
                call getCsr(mat, i, j, temp1)
                temp2 = temp2 + temp1 * vec(j)
            end do

            prod(i) = temp2

        end do
    end subroutine prodCsrVec



    ! ---------------------------------------------------------------------
    ! résolution Ax = b gauss seidel method
    ! ---------------------------------------------------------------------
    ! n : dimension de l'espace
    ! A : matrice creuse
    ! b : second membre
    ! eps : marge d'erreur (critère d'arret)
    ! x : vecteur solution retourné
    subroutine gaussSeidelCsr(n, A, b, eps, x)
        implicit none

        ! paramètres
        integer, intent(in) :: n
        type(CsrMat), intent(in) :: A
        real(rp), dimension(n), intent(in) :: b
        real(rp), intent(in) :: eps
        real(rp), dimension(n), intent(out) :: x

        ! variables locales
        integer :: i, k, cpt
        real(rp), dimension(n) :: xold
        real(rp) :: temp1, temp2, elt

        x = b
        ! pour rentrer dans le while
        xold = b + 1

        cpt = 0
        ! calcul de x itéré suivant
        do while (maxval(abs(x - xold)) > eps)
            cpt = cpt + 1
            xold = x

            do i = 1, n
                !call prodScal(i - 1, -A(i, 1:i-1), x(1:i-1), temp1)
                temp1 = 0.d0
                do k = A%jvcell(A%jposi(i)), A%jvcell(A%jposi(i + 1) - 1)
                    if (k > i - 1) then
                        exit
                    end if
                    call getCsr(A, i, k, elt)
                    temp1 = temp1 - elt * xold(k)
                end do

                !call prodScal(n - i, -A(i, i+1:n), xold(i+1:n), temp2)
                temp2 = 0.d0
                do k = A%jvcell(A%jposi(i)), A%jvcell(A%jposi(i + 1) - 1)
                    if (k > i) then
                        call getCsr(A, i, k, elt)
                        temp2 = temp2 - elt * xold(k)
                    end if
                end do


                call getCsr(A, i, i, elt)
                x(i) = (temp1 + temp2 + b(i)) / elt
            end do
        end do
        print *, "nb itérations gs    :", cpt
    end subroutine gaussSeidelCsr



    ! -------------------------------------------------------------------------------------------------------
    ! -------------------------------------------------------------------------------------------------------
    subroutine gradientConjugueCsr(A, b, eps, x)
        implicit none

        ! paramètres
        type(CsrMat), intent(in) :: A
        real(rp), dimension(A%n), intent(in) :: b
        real(rp), intent(in) :: eps
        real(rp), dimension(A%n), intent(out) :: x

        ! variables locales
        real(rp), dimension(A%n) :: r, rSuiv, xSuiv, p, pSuiv
        real(rp) :: alpha, beta
        real(rp) :: t1, t2
        real(rp), dimension(A%n) :: v1

        x = 0.0_rp
        call prodCsrVec(A, x, v1)
        r = b - v1
        p = r

        do while (maxval(abs(r)) > eps)
            call prodScal(A%n, r, r, t1)
            call prodCsrVec(A, p, v1)
            call prodScal(A%n, v1, p, t2)
            alpha = t1 / t2

            x = x + alpha * p

            r = r - alpha * v1

            call prodScal(A%n, r, r, t2)
            beta = t2 / t1

            p = r + beta * p
        end do
    end subroutine


    ! stocke la norme de X de taille N dans res
    !necessaire pour la condition d'arret de gauss seidel
    subroutine normeII(X,N,res)
        integer :: N,i
        double precision, dimension(N)::X
        double precision :: res
        res=0
        do i=1,N,1
            res=res+(X(i)**2)
        end do
        res=sqrt(res)
    end subroutine


    !resolution d'un systeme Ax=b par la methode de Gauss seidel
    ! recodé pour recevoir des matrice en CSR
    subroutine Gauss_Seidel_CSR(A,b,eps,resultat)
        type(CsrMat) :: A
        double precision, dimension(A%n) :: b, resultat,xk,xk1,test
        double precision :: eps,K1,K2,aii,test2
        integer :: i,j

        xk=1

        call prodCsrVec(A,xk,test)
        call normeII(test-b,A%n,test2)
        K1=0
        K2=0

        do while(test2>eps)
            do i=1,A%n,1
                do j=A%jposi(i), A%jposi(i+1)-1,1
                    if(A%jvcell(j)<i) then
                        K1=K1+A%tmat(j)*xk1(A%jvcell(j))
                    end if
                    if (A%jvcell(j)==i) then
                        aii=A%tmat(j)

                    end if
                    if (A%jvcell(j)>i) then
                        K2=K2+A%tmat(j)*xk(A%jvcell(j))
                    end if

                end do
                xk1(i)=(b(i)-K1-K2)/aii
                K1=0
                K2=0
            end do
            xk=xk1
            call prodCsrVec(A,xk,test)
            call normeII(test-b,A%n,test2)
        end do
        resultat = xk
    end subroutine


END MODULE csr
