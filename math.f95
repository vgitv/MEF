! =====================================================================
! Module contenant les fonctions :
!
! linspace
! seq
! evaluate
! evaluate2D
! saveSol
! saveMat
! affMat
! id
! plu
! cholesky
! prodScal
! prodMatVec
! descenteRemonte
! linSolve
! norme2
! normeInf
! trapeze
! =====================================================================

MODULE math

    implicit none

    integer, parameter :: sp = kind(1.0E0)
    integer, parameter :: dp = kind(1.0D0)
    integer, parameter :: rp = dp
    real(rp), parameter :: pi = 4.0_rp * atan(1.0_rp)

contains

    ! ---------------------------------------------------------------------
    ! créer un vecteur uniforme en précisant la taille
    ! ---------------------------------------------------------------------
    ! a : valeur de début
    ! b : valeur de fin
    ! n : nombre de valeur du vecteur retourné
    function linspace(a, b, n)
        implicit none

        ! param
        real(rp), intent(in) :: a, b
        integer, intent(in) :: n

        ! return
        real(rp), dimension(n) :: linspace

        ! variables locales
        real(rp) :: h
        integer :: i

        h = (b - a) / (n - 1)
        !linspace(1:n) = (a + i * h, i = 0, n - 1)
        do i = 0, n - 2
            linspace(i + 1) = a + i * h
        end do
        linspace(n) = b
    end function linspace



    ! ---------------------------------------------------------------------
    ! créer un vecteur uniforme en précisant le pas
    ! ---------------------------------------------------------------------
    ! a : valeur de début
    ! b : valeur de fin
    ! h : pas éventuellement diminué pour que la valeur finale soit b
    function seq(a, b, h)
        implicit none

        ! param
        real(rp), intent(in) :: a, b
        real(rp) :: h

        ! return
        real(rp), dimension(:), allocatable :: seq

        ! variables locales
        integer :: n, i

        ! nombre de valeurs
        n = ceiling((b - a) / h) + 1
        allocate(seq(n))

        ! nouveau pas inf ou égal au param
        h = (b - a) / (n - 1)
        seq = linspace(a, b, n)
    end function seq



    ! ---------------------------------------------------------------------
    ! evalue une fonction sur un vecteur
    ! ---------------------------------------------------------------------
    ! fc : fonction de une variable
    ! vec : vecteur quelconque
    ! n : taille de vec
    !
    ! sortie : fc(vec) qui est un vecteur de taille n
    function evaluate(fc, vec)
        implicit none

        ! param
        real(rp), external :: fc
        real(rp), dimension(:), intent(in) :: vec

        ! return
        real(rp), dimension(size(vec)) :: evaluate

        ! variables locales
        integer :: i

        do i = 1, size(vec)
            evaluate(i) = fc(vec(i))
        end do
    end function

    function evaluate2D(fc, vx, vy)
        implicit none

        ! PARAMÈTRES
        real(rp), external :: fc
        real(rp), dimension(:), intent(in) :: vx
        real(rp), dimension(:), intent(in) :: vy

        ! RETURN
        real(rp), dimension(size(vx)) :: evaluate2D

        ! VARIABLES LOCALES
        integer :: i

        do i = 1, size(vx)
            evaluate2D(i) = fc(vx(i), vy(i))
        end do
    end function



    ! ---------------------------------------------------------------------
    ! sauvegarde les antécédants et la solution en colonne dans un fichier
    ! ---------------------------------------------------------------------
    ! x : abscisse
    ! sol : ordonnée
    ! n : taille de x et sol
    ! fichier : nom du fichier créé
    subroutine saveSol(x, sol, fichier)
        implicit none

        ! parametres
        real(rp), dimension(:) :: x, sol
        character(len=*) :: fichier

        ! variables locales
        integer :: i

        open(unit=0, file=fichier)
        do i = 1, size(x)
            write (0, *) x(i), sol(i)
        end do
        close(0)
    end subroutine saveSol



    ! ---------------------------------------------------------------------
    ! sauvegarde les antécédants et la solution en colonne dans un fichier
    ! ---------------------------------------------------------------------
    ! x : abscisse
    ! mat : matrice de n lignes
    ! n : taille de x
    ! fichier : nom du fichier créé
    subroutine saveMat(x, mat, n, fichier)
        implicit none

        ! parametres
        real(rp), dimension(:) :: x
        real(rp), dimension(:, :) :: mat
        integer :: n
        character(len=*) :: fichier

        ! variables locales
        integer :: i

        open(unit=0, file=fichier)
        do i = 1, n
            write (0, *) x(i), mat(i, :)
        end do
        close(0)
    end subroutine saveMat



    ! ---------------------------------------------------------------------
    ! affichage matrice en console
    ! ---------------------------------------------------------------------
    ! n : nombre de ligne de la matrice
    ! mat : matrice à afficher
    subroutine affMat(mat)
        implicit none

        ! parametres
        real(rp), dimension(:, :), intent(in) :: mat

        ! variables locales
        integer :: i, n
        integer, dimension(2) :: forme

        forme = shape(mat)
        n = forme(1)

        do i = 1, n
            print "(100F8.2)", mat(i, :)
        end do
    end subroutine affMat



    ! ---------------------------------------------------------------------
    ! matrice identité
    ! ---------------------------------------------------------------------
    ! n : taille de la matrice à créer
    ! mat : matrice identité n par n retournée
    subroutine id(n, mat)
        implicit none

        ! parametres
        integer, intent(in) :: n
        real(rp), dimension(n, n), intent(out) :: mat

        ! variables locales
        integer :: i

        mat = 0.d0
        do i = 1, n
            mat(i, i) = 1.d0
        end do
    end subroutine id



    ! ---------------------------------------------------------------------
    ! décomposition PLU
    ! ---------------------------------------------------------------------
    ! n : dimension de l'espace
    ! A : matrice à décomposer
    !
    ! retours :
    ! P : matrice d'interversion ligne si pivot nul
    ! L : triang inf
    ! U : triang sup
    subroutine plu(n, A, P, L, U)
        implicit none

        ! parametres
        integer, intent(in) :: n
        real(rp), dimension(n, n), intent(in) :: A
        real(rp), dimension(n, n), intent(out) :: P, L, U

        ! variables locales
        integer :: i, j, k
        real(rp), dimension(n) :: temp
        real(rp) :: maxk

        call id(n, P)
        call id(n, L)
        U = A

        do k = 1, n - 1

            ! détermination de i tq A(i, k) = max_j(A(j, k)) j>=k
            i = k
            maxk = abs(A(i, k))
            do j = k + 1, n
                if (abs(A(j, k)) > maxk) then
                    maxk = abs(A(j, k))
                    i = j
                end if
            end do

            ! permutation lignes colonnes
            if (i /= k) then
                ! échange des lignes k et i de P, L, U
                P((/k, i/), :) = P((/i, k/), :)
                L((/k, i/), :) = L((/i, k/), :)
                U((/k, i/), :) = U((/i, k/), :)

                ! échange des colonnes k et i de L
                L(:, (/k, i/)) = L(:, (/i, k/))
            end if

            do i = k + 1, n
                L(i, k) = U(i, k) / U(k, k)
                do j = k + 1, n
                    U(i, j) = U(i, j) - L(i, k) * U(k, j)
                end do
                U(i, k) = 0.d0
            end do

        end do
    end subroutine plu



    ! ---------------------------------------------------------------------
    ! décomposition Cholesky (cas A SDP)
    ! ---------------------------------------------------------------------
    ! n : dimension d'espace
    ! A : matrice à décomposer
    ! L : matrice symétrique issue de la décomposition (n * n)
    subroutine cholesky(n, A, L)
        implicit none

        ! paramètres
        integer, intent(in) :: n
        real(rp), dimension(n, n), intent(in) :: A
        real(rp), dimension(n, n), intent(out) :: L

        ! variables locales
        integer :: i, j

        L(1, 1) = sqrt(A(1, 1))

        do i = 2, n
            L(i, 1) = A(i, 1) / L(1, 1)
            L(i, i) = sqrt(A(i, i) - sum(L(i, 1:i-1)**2))
            do j = i + 1, n
                L(j, i) = (A(j, i) - sum(L(i, 1:i-1) * L(j, 1:i-1))) / L(i, i)
            end do
        end do
    end subroutine cholesky



    ! ---------------------------------------------------------------------
    ! produit vecteur vecteur
    ! ---------------------------------------------------------------------
    ! n : dimension de l'espace
    ! vec1, vec2 : deux vecteurs de taille n
    ! prod : scalaire, résultat du produit
    subroutine prodScal(n, vec1, vec2, prod)
        implicit none

        ! parametres
        integer, intent(in) :: n
        real(rp), dimension(n), intent(in) :: vec1, vec2
        real(rp), intent(out) :: prod

        ! variables locales
        integer :: i

        prod = 0.d0

        do i = 1, n
            prod = prod + vec1(i) * vec2(i) 
        end do
    end subroutine prodScal



    ! ---------------------------------------------------------------------
    ! produit matrice vecteur
    ! ---------------------------------------------------------------------
    ! n : dimension de l'espace
    ! mat : matrice n*n
    ! vec : vecteur de taille n
    ! prod : vecteur de taille n résultat du produit
    subroutine prodMatVec(n, mat, vec, prod)
        implicit none

        ! parametres
        integer, intent(in) :: n
        real(rp), dimension(n, n), intent(in) :: mat
        real(rp), dimension(n), intent(in) :: vec
        real(rp), dimension(n), intent(out) :: prod

        ! variables locales
        integer :: i
        real(rp) :: temp

        do i = 1, n
            call prodScal(n, mat(i, :), vec, temp)
            prod(i) = temp
        end do
    end subroutine prodMatVec



    ! ---------------------------------------------------------------------
    ! résolution x tq Ly = Pb avec Ux = y
    ! ---------------------------------------------------------------------
    ! n : dimension de l'espace
    ! b : second membre
    ! P : passage
    ! L : triang ing
    ! U : triang sup
    ! x : solution retournée
    subroutine descenteRemontee(n, b, P, L, U, x)
        implicit none

        ! parametres
        integer, intent(in) :: n! dimension de l'espace
        real(rp), dimension(n), intent(in) :: b! second menbre
        real(rp), dimension(n, n), intent(in) :: P, L, U! issu décomposition plu
        real(rp), dimension(n), intent(out) :: x! solution

        ! variables lodales
        integer :: i
        real(rp), dimension(n) :: Pb, verif
        real(rp) :: temp
        real(rp), dimension(n) :: y

        call prodMatVec(n, P, b, Pb)

        ! algorithme de desente y tq Ly = Pb
        y(1) = Pb(1) / L(1, 1)

        do i = 2, n
            call prodScal(i - 1, y(1:i-1), L(i, 1:i-1), temp)
            y(i) = (Pb(i) - temp) / L(i, i)
        end do

        ! algorithme de remontée x tq Ux = y
        print *, U(n, n), y(n)
        x(n) = y(n) / U(n, n)

        do i = n - 1, 1, -1
            call prodScal(n - i, x(i+1:n), U(i, i+1:n), temp)
            x(i) = (y(i) - temp) / U(i, i)
        end do

        call prodMatVec(n, L, y, verif)
    end subroutine descenteRemontee



    ! ---------------------------------------------------------------------
    ! résolution Ax = b
    ! ---------------------------------------------------------------------
    ! method : au choix parmi plu, cholesky
    ! n : dimension de l'espace
    ! A : matrice n*n
    ! b : second membre
    ! x : solution retournée
    subroutine linSolve(method, n, A, b, x)
        implicit none

        ! parametres
        character(len = *), intent(in) :: method! plu ou cholesky
        integer, intent(in) :: n! dimension de l'espace
        real(rp), dimension(n, n), intent(in) :: A! application linéaire
        real(rp), dimension(n), intent(in) :: b! second membre
        real(rp), dimension(n), intent(out) :: x! solution

        ! variables locales
        real(rp), dimension(n, n) :: P, L, U

        if (method == "plu") then
            call plu(n, A, P, L, U)
        else if (method == "cholesky") then
            call id(n, P)
            call cholesky(n, A, L)
            U = transpose(L)
        else
            print *, "choix de méthode invalide"
        end if

        call descenteRemontee(n, b, P, L, U, x)
    end subroutine linSolve



    ! ---------------------------------------------------------------------
    ! norme euclidienne discrète
    ! ---------------------------------------------------------------------
    ! v : vecteur
    ! norme : valeur retournée
    subroutine norme2(v, norme)
        implicit none

        ! paramètres
        real(rp), dimension(:), intent(in) :: v
        real(rp), intent(out) :: norme

        norme = sqrt(sum(v**2))
    end subroutine norme2



    ! ---------------------------------------------------------------------
    ! norme infinie discrète
    ! ---------------------------------------------------------------------
    ! v : vecteur
    ! norme : valeur retournée
    subroutine normeInf(v, norme)
        implicit none

        ! paramètres
        real(rp), dimension(:), intent(in) :: v
        real(rp), intent(out) :: norme

        norme = maxval(abs(v))
    end subroutine normeInf



    ! ---------------------------------------------------------------------
    ! méthode des trapèzes
    ! ---------------------------------------------------------------------
    ! x : vecteur des antécédants (subdivision abscisse)
    ! fx : vecteur des images (fx = f(x))
    ! somme : approximation de l'intégrale retournée
    subroutine trapeze(x, fx, somme)
        implicit none

        ! paramètres
        real(rp), dimension(:), intent(in) :: x, fx
        real(rp), intent(out) :: somme

        ! variables locales
        integer :: i, n

        somme = 0.d0
        n = size(x)

        do i = 1, n - 1
            somme = somme + (x(i + 1) - x(i)) * (fx(i) + fx(i + 1)) / 2
        end do
    end subroutine trapeze



    subroutine saveRMat(mat, fichier)
        implicit none

        ! paramètres
        real(rp), dimension(:, :) :: mat
        character(len = *), intent(in) :: fichier

        ! variables locales
        integer, dimension(2) :: forme
        integer :: i

        forme = shape(mat)

        open(unit = 1, file = fichier)
        do i = 1, forme(1)
            write(1, *), mat(i, :)
        end do
        close(1)
    end subroutine

    subroutine saveIMat(mat, fichier)
        implicit none

        ! paramètres
        integer, dimension(:, :) :: mat
        character(len = *), intent(in) :: fichier

        ! variables locales
        integer, dimension(2) :: forme
        integer :: i

        forme = shape(mat)

        open(unit = 1, file = fichier)
        do i = 1, forme(1)
            write(1, *), mat(i, :)
        end do
        close(1)
    end subroutine

    subroutine trieBulle(tab)
        implicit none

        ! paramètres
        integer, dimension(:) :: tab

        ! variables locales
        integer :: i, j, temp, imax
        integer :: n

        n = size(tab)

        do j = n, 1, -1
            imax = 1
            temp = tab(1)
            do i = 1, j
                ! si on trouve un elt plus grand
                if (tab(i) > temp) then
                    imax = i
                    ! c'est le nouveau plus grand elt
                    temp = tab(i)
                end if
            end do

            ! on inverse le plus grand elt et le dernier elt
            temp = tab(imax)
            tab(imax) = tab(j)
            tab(j) = temp
        end do
    end subroutine

END MODULE math
