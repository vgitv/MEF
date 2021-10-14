! ===========================================================================================================
! EF P1 lagrange
!
! équation de poisson 2D
! condition de dirichlet non homogène
!
! Fichiers sources :
! main.f95
! elliptique.f95
! csr.f95
! math.f95
! sdd.f95
! donnees.f95
! ===========================================================================================================

PROGRAM main

    use donnees
    use math
    use sdd
    use csr
    use elliptique

    implicit none

    type(Maillage) :: m
    real(rp), dimension(:), allocatable :: u, uOmega, sol
    type(CsrMat) :: A
    real(rp), dimension(:), allocatable :: l
    real(rp) :: eps = 0.001
    integer :: i
    logical :: aff = .false.
    real(rp) :: t1, t2



    ! =======================================================================================================
    ! partie principale
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! construction maillage
    ! -------------------------------------------------------------------------------------------------------
    ! lecture maillage fichier .mesh
    call readMesh("maillage/out.final.mesh", m)

    ! liste des segments (par listes chaînées)
    call findSeg(m)

    ! écriture du maillage fichiers pour vérif
    call exportMesh(m, "sorties")



    ! -------------------------------------------------------------------------------------------------------
    ! construction matrice de masse globale
    ! -------------------------------------------------------------------------------------------------------
    ! assemblage matrice de masse globale
    call buildGlobalMatrix(A, m, d_phi1, d_phi2, d_phi3)
    if (aff) then
        ! affichage pour vérifications, désactivés par défaut
        print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print *, "jposi ", A%jposi
        print *
        print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print *, "jvcell", A%jvcell
        print *
        print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print *, "tmat", A%tmat
        print *
        print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print *, "matrice de masse globale"
        call affCsr(A)
    end if



    ! -------------------------------------------------------------------------------------------------------
    ! assemblage second membre
    ! -------------------------------------------------------------------------------------------------------
    call buildL(m, f, t0, L)
    if (aff) then
        ! affichage pour vérifications, désactivés par défaut
        print *
        print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print *, "second membre L", L
    end if



    ! -------------------------------------------------------------------------------------------------------
    ! résolution système
    ! -------------------------------------------------------------------------------------------------------
    allocate(u(m%ni))


    ! GAUSS SEIDEL PAS EFFICACE
    !call cpu_time(t1)
    !call gaussSeidelCsr(A%n, A, L, eps, u)
    !call cpu_time(t2)
    !print *, "gauss seidel time :", t2 - t1

    call cpu_time(t1)
    call gradientConjugueCsr(A, L, eps, u)
    call cpu_time(t2)
    print *, "gradient conjugue time :", t2 - t1



    ! -------------------------------------------------------------------------------------------------------
    ! affectation valeurs au bord
    ! -------------------------------------------------------------------------------------------------------
    allocate(uOmega(m%ns))
    call edge(m, t0, u, uOmega)
    if (aff) then
        ! affichage pour vérifications, désactivés par défaut
        print *
        print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print *, "sol approchée uOmega", uOmega
    end if
    call CellVertexVtk(m, uOmega, "sorties/solution.vtk")



    ! -------------------------------------------------------------------------------------------------------
    ! construction solution exacte pour cas vérifiables
    ! -------------------------------------------------------------------------------------------------------
    ! solution exacte
    allocate(sol(m%ns))
    call buildSol(m, uExacte, sol)
    call CellVertexVtk(m, sol, "sorties/test.vtk")
    if (aff) then
        ! affichage pour vérifications, désactivés par défaut
        print *
        print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print *, "solution exacte     ", sol
    end if



    ! -------------------------------------------------------------------------------------------------------
    ! messages console
    ! -------------------------------------------------------------------------------------------------------
    print *
    print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print *, "Erreur relative norme infinie", maxval(abs(uOmega - sol)) / maxval(abs(sol))
    print *, "Erreur relative norme L2     ", sqrt(sum((uOmega - sol)**2)) / sqrt(sum(uOmega))
    print *, "max", maxval(uOmega)
    print *, "min", minval(uOmega)






    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DÉSALLOCATIONS FINALES
    call rmMaillage(m)
    call rmCsr(A)
    deallocate(L, u, sol)
END PROGRAM main
