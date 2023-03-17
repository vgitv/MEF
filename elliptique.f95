MODULE elliptique

    use math
    use sdd
    use csr

    implicit none

    type :: Maillage
        ! nombre de sommets
        integer :: ns, ni, nb
        ! m%coordonnées des sommets
        real(rp), dimension(:, :), allocatable :: coord
        ! est-ce un terme de bord ? + numérotation locale des sommets intérieurs
        integer, dimension(:), allocatable :: bord


        ! nombre d'éléments
        integer :: ne
        ! numérotation globale de chaque sommet
        integer, dimension(:, :), allocatable :: connect


        ! nombre de sommets
        integer :: nseg
        ! liste des segments
        integer, dimension(:, :), allocatable :: nubo
    end type Maillage

contains

    ! =======================================================================================================
    ! FONCTION RELATIVES AU MAILLAGE
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! destructeur Maillage
    ! -------------------------------------------------------------------------------------------------------
    subroutine rmMaillage(m)
        implicit none

        type(Maillage) :: m

        if (allocated(m%coord)) then
            deallocate(m%coord)
        end if

        if (allocated(m%bord)) then
            deallocate(m%bord)
        end if

        if (allocated(m%connect)) then
            deallocate(m%connect)
        end if

        if (allocated(m%nubo)) then
            deallocate(m%nubo)
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! constructeur Maillage : lit le fichier .mesh et renvoie toutes les infos sur le maillage
    ! -------------------------------------------------------------------------------------------------------
    ! fichier : le chemin vers le fichier *.mesh
    !
    ! SORTIE :
    ! ns       : nombre de sommets
    ! m%coord    : tableau des m%coordonnées cartésiennes des sommet
    ! bord     : est-ce un terme de bord ? (0 / 1)
    ! ne       : nombre d'éléments
    ! noeudGeo : numérotation globale des sommets
    subroutine readMesh(fichier, m)!ns, m%coord, bord, ne, connect)
        implicit none

        ! PARAMÈTRES
        character(len = *), intent(in) :: fichier
        type(Maillage), intent(inout) :: m

        ! VARIABLES LOCALES
        integer :: i, cpt, loc, last

        open(unit = 0, file = fichier)
        ! on saute les 4 premières lignes
        do i = 1, 4
            read (0, *)
        end do

        ! lecture nombre de sommets
        read (0, *) m%ns
        allocate(m%coord(m%ns, 2), m%bord(m%ns))

        ! lecture m%coordonnées des sommets et si ils sont sur le bord
        loc = -1
        last = 1
        m%ni = 0
        do cpt = 1, m%ns
            read (0, *) m%coord(cpt, :), m%bord(cpt)
            if (m%bord(cpt) == 0) then
                m%ni = m%ni + 1
                m%bord(cpt) = loc
                loc = loc - 1
            else
                m%bord(cpt) = last
                last = last + 1
            end if
        end do

        ! on saute une ligne
        read (0, *)
        ! lecture du nombre de m%connects
        read (0, *) m%ne
        allocate(m%connect(m%ne, 3))

        ! lecture numérotation des sommets
        do cpt = 1, m%ne
            read (0, *) m%connect(cpt, :)
        end do
        close(0)

        m%nb = m%ns - m%ni
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! créé la liste des segment nubo
    ! -------------------------------------------------------------------------------------------------------
    subroutine findSeg(mesh)
        implicit none

        ! paramètres
        type(Maillage), intent(inout) :: mesh

        ! variables locales
        type(Liste), dimension(mesh%ns) :: hashtab
        integer :: i, k, s
        integer :: ismin, ismax
        integer :: suiv
        integer :: nseg
        logical :: insert
        type(Cell), pointer :: p, it

        ! initialisation de hashtab
        do i = 1, mesh%ns
            call newListe(hashtab(i))
            call insertListe(hashtab(i), NULL(), i)
        end do

        ! initialisation de nseg
        nseg = 0

        ! construction de hashtab

        ! parcourt des triangles K
        do k = 1, mesh%ne

            ! parcourt des sommets du triangles
            do s = 1, 3
                ! récupération indices globaux
                suiv = s + 1
                if (suiv > 3) then
                    suiv = 1
                end if

                if (mesh%connect(k, s) < mesh%connect(k, suiv)) then
                    ismin = mesh%connect(k, s)
                    ismax = mesh%connect(k, suiv)
                else
                    ismin = mesh%connect(k, suiv)
                    ismax = mesh%connect(k, s)
                end if

                ! insertion
                insert = insertTrieListe(hashtab(ismin), ismax)
                if (insert) then
                    nseg = nseg + 1
                end if
            end do
        end do

        ! construction maillage
        mesh%nseg = nseg
        k = 1
        allocate(mesh%nubo(nseg, 2))
        do i = 1, mesh%ns
            p => hashtab(i)%premier
            if (associated(p%suiv)) then
                it => p%suiv
                do while (associated(it))
                    mesh%nubo(k, 1) = p%elt
                    mesh%nubo(k, 2) = it%elt
                    it => it%suiv
                    k = k + 1
                end do
            end if
        end do

        ! destruction de hashtab
        do i = 1, mesh%ns
            call rmListe(hashtab(i))
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! évalu les fonctions de base aux points de la quadrature de gauss
    ! -------------------------------------------------------------------------------------------------------
    subroutine evalBase(mesh, e1, e2, e3)
        implicit none

        ! paramètres
        type(Maillage), intent(inout) :: mesh
        real(rp), external :: e1, e2, e3

        ! variables locales
        !++!
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! sort un fichier au format .vtk
    ! -------------------------------------------------------------------------------------------------------
    ! m : maillage
    ! solution : solution du problème sur les points du maillage
    subroutine CellVertexVtk(m, solution, fichier)

        type(Maillage), intent(in) :: m
        real(rp), dimension(m%ns), intent(in) :: solution
        character(len = *), intent(in) :: fichier

        integer :: lpbName, lensd3, i
        character(len=75) :: vtk



        !lpbName = INDEX(pbName,' ')-1
        !vtk(1:lpbName) = pbName(1:lpbName)
        !vtk(lpbName+1:lpbName+5)='.vtk'
        !lensd3 =lpbName+4
        OPEN(unit = 10, file = fichier, action = 'write', status = 'replace')
        write(10, '(A)') '# vtk DataFile Version 3.0'
        write(10, '(A)') '# solution du M1 TR'
        write(10, '(A)') 'ASCII'
        write(10, '(A)') 'DATASET UNSTRUCTURED_GRID'
        write(10, '(A,I7,A)') 'POINTS ', m%ns,' float'
        !----------------------------------------------------------------------
        Do i = 1, m%ns
            !write(10,'(E13.7,2x,E13.7,2x,E13.7)') m%coord(i, 1), m%coord(i, 2), 0.0_rp
            write(10,'(E13.7,2x,E13.7,2x,E13.7)') m%coord(i, 1), m%coord(i, 2), solution(i)
        end do
        !-----------------------------------------------------------------------
        write(10, '(A,1x,I7,1x,I6)')'CELLS',m%Ne,4 * m%Ne

        Do i = 1, m%Ne
            write(10,'(I7,1x,I7,1x,I7)') 3, m%connect(i, 1) - 1, m%connect(i, 2) - 1, m%connect(i, 3)-1
        end do

        write(10, '(A,1x,I7)') 'CELL_TYPES', m%Ne

        Do i = 1, m%Ne
            write(10, '(I1)') 5
        end do
        !---------------------------------------------------------------------
        write(10, '(A,1x,I7)') 'POINT_DATA', m%ns
        write(10, '(A)') 'SCALARS hauteur float'
        write(10, '(A)') 'LOOKUP_TABLE DEFAULT'

        Do i = 1, m%ns
            write(10,'(ES20.7)') solution(i)
        end do

        !-------------------------------------------------------------------
        write(10,'(A)') 'SCALARS topographie float'
        write(10, '(A)') 'LOOKUP_TABLE DEFAULT'

        Do i = 1, m%ns
            write(10,'(ES20.7)') max(1.E-08,0.0)
        end do
        !-----------------------------------------------------------------------
        write(10,'(A)') 'VECTORS vitesse float'

        Do i = 1, m%ns
            write(10, '(E13.7,2x,E13.7,2x,E13.7)') 1.E-08, 1.E-08,1.E-08
        end do

        close(10)
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! sauvegarde les objet dans un format lisible par Gnuplot
    ! -------------------------------------------------------------------------------------------------------
    subroutine exportMesh(m, dossier)
        implicit none

        ! paramètres
        type(Maillage), intent(in) :: m
        character(len = *), intent(in) :: dossier

        ! variables locales
        integer :: i

        call saveRMat(m%coord, dossier // "/coord.dat")

        open(unit = 1, file = dossier // "/bord.dat")
        do i = 1, m%ns
            write(1, *), m%bord(i)
        end do
        close(1)

        call saveIMat(m%connect, dossier // "/connect.dat")
        call saveIMat(m%nubo, dossier // "/nubo.dat")

        open(unit = 1, file = dossier // "/n.dat")
        write(1, *), "ns   :", m%ns
        write(1, *), "ni   :", m%ni
        write(1, *), "nb   :", m%nb
        write(1, *), "ne   :", m%ne
        write(1, *), "nseg :", m%nseg
        close(1)
    end subroutine



    ! =======================================================================================================
    ! FONCTIONS RELATIVES À L'ASSENBLAGE DE A
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! calcule le nombre de valeur non nulle de la matrice de masse globale, à partir de la liste des segments
    ! -------------------------------------------------------------------------------------------------------
    ! A : future matrice de masse globale
    ! mesh : maillage
    subroutine nCoeffMat(A, mesh)
        implicit none

        ! paramètres
        type(CsrMat), intent(inout) :: A
        type(Maillage), intent(in) :: mesh

        ! variables locales
        integer :: nCoeff
        integer :: iseg, is, js

        nCoeff = mesh%ni

        do iseg = 1, mesh%nseg
            is = mesh%nubo(iseg, 1)
            js = mesh%nubo(iseg, 2)

            ! si on est sur un segment intérieur
            if ((mesh%bord(is) < 0) .and. (mesh%bord(js) < 0)) then
                nCoeff = nCoeff + 2
            end if
        end do

        A%nCoeff = nCoeff
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! calcule jPosi, vecteur qui compte les 1ers elts de chaque ligne d'une matrice csr
    ! -------------------------------------------------------------------------------------------------------
    ! A : future matrice de masse globale
    ! mesh : maillage
    subroutine jposiMat(A, mesh)
        implicit none

        ! paramètres
        type(CsrMat), intent(inout) :: A
        type(Maillage), intent(in) :: mesh

        ! variables locales
        integer :: i
        integer :: is, js, iseg
        integer :: isLocal, jsLocal

        ! éléments diagonaux
        do i = 1, A%n + 1
            A%jposi(i) = i
        end do

        ! on se sert de cet indice inutilisé pour stocker la taille de jposi (évite de recalculer ns + 1)
        i = A%n + 1
        do iseg = 1, mesh%nseg
            is = mesh%nubo(iseg, 1)
            js = mesh%nubo(iseg, 2)

            ! si on est sur un segment intérieur
            if ((mesh%bord(is) < 0) .and. (mesh%bord(js) < 0)) then
                !isLocal = is - sum(mesh%bord(1 : is))
                !jsLocal = js - sum(mesh%bord(1 : js))
                isLocal = -mesh%bord(is)
                jsLocal = -mesh%bord(js)

                A%jposi(isLocal + 1 : i) = A%jposi(isLocal + 1 : i) + 1
                A%jposi(jsLocal + 1 : i) = A%jposi(jsLocal + 1 : i) + 1
            end if
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! calcule jvcell, vecteur qui compte les 1ers elts de chaque ligne d'une matrice csr
    ! -------------------------------------------------------------------------------------------------------
    ! A : future matrice de masse globale
    ! mesh : maillage
    subroutine jvcellMat(A, mesh)
        implicit none

        ! paramètres
        type(CsrMat), intent(inout) :: A
        type(Maillage), intent(in) :: mesh

        ! variables locales
        integer, dimension(A%n) :: indPL
        integer :: i, iseg, is, js
        integer :: isLocal, jsLocal

        ! initialisation de indPL
        indPL(1 : A%n) = A%jposi(1 : A%n)

        ! indices de colonne des termes diagonaux
        do i = 1, A%n
            A%jvcell(indPL(i)) = i
            indPL(i) = indPL(i) + 1
        end do

        ! autres indices
        do iseg = 1, mesh%nseg
            is = mesh%nubo(iseg, 1)
            js = mesh%nubo(iseg, 2)

            ! si on est sur un segment intérieur
            if ((mesh%bord(is) < 0) .and. (mesh%bord(js) < 0)) then
                !isLocal = is - sum(mesh%bord(1 : is))
                !jsLocal = js - sum(mesh%bord(1 : js))
                isLocal = -mesh%bord(is)
                jsLocal = -mesh%bord(js)

                A%jvcell(indPL(isLocal)) = jsLocal
                indPL(isLocal) = indPL(isLocal) + 1

                A%jvcell(indPL(jsLocal)) = isLocal
                indPL(jsLocal) = indPL(jsLocal) + 1
            end if
        end do

        ! on trie chacun de nos sous tableaux
        ! le premier sous tableau est déjà trié
        !call trieBulle(A%jvcell(1 : indPL(1) - 1))
        do i = 1, A%n - 1
            call trieBulle(A%jvcell(indPL(i) : indPL(i + 1) - 1))
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! construction contribution
    ! -------------------------------------------------------------------------------------------------------
    ! k : numéro du triangle (élément)
    ! mesh : maillage
    ! d1, d2, d3 : gradient des fonctions de base
    function AmatLoc(k, mesh, d1, d2, d3)
        implicit none

        ! paramètres
        integer, intent(in) :: k
        type(Maillage), intent(in) :: mesh
        real(rp), dimension(2) :: d1, d2, d3

        ! return
        real(rp), dimension(3, 3) :: AmatLoc

        ! variables locales
        real(rp) :: x1, y1, x2, y2, x3, y3
        real(rp), dimension(2, 2) :: Jmoins1
        real(rp) :: detJ
        real(rp), dimension(3, 2) :: v
        integer :: i, j
        real(rp) :: a, b, c, p, aire, temp

        ! coordonnées du premier point
        x1 = mesh%coord(mesh%connect(k, 1), 1)
        y1 = mesh%coord(mesh%connect(k, 1), 2)

        ! coordonnées du deuxième point
        x2 = mesh%coord(mesh%connect(k, 2), 1)
        y2 = mesh%coord(mesh%connect(k, 2), 2)

        ! coordonnées du troisième point
        x3 = mesh%coord(mesh%connect(k, 3), 1)
        y3 = mesh%coord(mesh%connect(k, 3), 2)

        ! aire du triangle formule de héron
        !call norme2((/x1, y1/) - (/x2, y2/), a)
        !call norme2((/x2, y2/) - (/x3, y3/), b)
        !call norme2((/x3, y3/) - (/x1, y1/), c)
        !p = (a + b + c) / 2
        !aire = sqrt(p * (p - a) * (p - b) * (p - c))
        ! aire du triangle de référence
        aire = 0.5_rp

        ! déterminant de la jacobienne transformation K ref -> K
        detJ = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)
        !print *, "déterminant J", detJ

        ! inverse de la jacobienne (formule des cofacteurs)
        Jmoins1(1, 1) = (y3 - y1) / detJ
        Jmoins1(1, 2) = (x1 - x3) / detJ
        Jmoins1(2, 1) = (y1 - y2) / detJ
        Jmoins1(2, 2) = (x2 - x1) / detJ

        ! les deux termes de l'intégrale
        call prodMatVec(2, Jmoins1, d1, v(1, :))
        call prodMatVec(2, Jmoins1, d2, v(2, :))
        call prodMatVec(2, Jmoins1, d3, v(3, :))

        ! construction matrice locale 3 par 3
        do i = 1, 3
            do j = 1, 3
                ! approximation de l'intégrale de vi.vj
                call prodScal(2, v(i, :), v(j, :), temp)
                AmatLoc(i, j) = temp * abs(detJ) * aire
            end do
        end do

        !print *, "### AmatLoc Triangle ", k, "sommets ", mesh%connect(k, :)
        !call affMat(AmatLoc)
        !print *
        !print *
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! ajout la contribution locale à A en stockage Csr
    ! -------------------------------------------------------------------------------------------------------
    subroutine ajout(ii, jj, coefA, A)
        implicit none

        ! paramètres
        integer, intent(in) :: ii
        integer, intent(in) :: jj
        real(rp), intent(in) :: coefA
        type(CsrMat), intent(inout) :: A

        ! variables locales
        logical :: trouver
        integer :: j

        trouver = .false.
        ! parcourt de la ligne ii
        do j = A%jposi(ii), A%jposi(ii + 1) - 1
            ! trouver le bon numéro de colonne
            if (A%jvcell(j) == jj) then
                A%tmat(j) = A%tmat(j) + coefA
                trouver = .true.
                exit
            end if
        end do

        if (.not. trouver) then
            print *, "problème d'assemblage de la matrice A"
            stop
        end if
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! ajout des 9 contributions dans la matrice de masse globale
    ! -------------------------------------------------------------------------------------------------------
    subroutine assemble(mesh, k, Aloc, A)
        implicit none

        ! paramètres
        type(Maillage), intent(in) :: mesh
        integer, intent(in) :: k
        real(rp), dimension(3, 3), intent(in) :: Aloc
        type(CsrMat), intent(inout) :: A

        ! variables locales
        integer, dimension(3) :: s
        integer :: i, j, ilocal, jlocal

        ! numéros globaux des sommets du triangle k
        s(1) = mesh%connect(k, 1)
        s(2) = mesh%connect(k, 2)
        s(3) = mesh%connect(k, 3)

        ! ligne i
        do i = 1, 3
            ! colonne j
            do j = 1, 3
                if ((mesh%bord(s(i)) < 0) .and. (mesh%bord(s(j)) < 0)) then
                    ilocal = -mesh%bord(s(i))
                    jlocal = -mesh%bord(s(j))
                    call ajout(ilocal, jlocal, Aloc(i, j), A)
                end if
            end do
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! calcule tmat, vecteur des valeurs non nulles de la matrice de masse globale
    ! -------------------------------------------------------------------------------------------------------
    ! A : future matrice de masse globale
    ! mesh : maillage
    ! d1, d2, d3 : gradient des fonctions de base (cst pour P1 lagrange)
    subroutine tmatMat(A, mesh, d1, d2, d3)
        implicit none

        ! paramètres
        type(CsrMat), intent(inout) :: A
        type(Maillage), intent(in) :: mesh
        real(rp), dimension(2), intent(in) :: d1, d2, d3

        ! variables globales
        integer :: k
        real(rp), dimension(3, 3) :: Aloc

        A%tmat = 0.0_rp

        do k = 1, mesh%ne
            Aloc = AmatLoc(k, mesh, d1, d2, d3)
            call assemble(mesh, k, Aloc, A)
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! assemblage de la matrice de masse globale
    ! -------------------------------------------------------------------------------------------------------
    ! A : future matrice de masse globale
    ! mesh : maillage
    subroutine buildGlobalMatrix(A, mesh, d1, d2, d3)
        implicit none

        ! paramètres
        type(CsrMat), intent(out) :: A
        type(Maillage), intent(in) :: mesh
        real(rp), dimension(2) :: d1, d2, d3

        ! matrice carré, n = m
        A%n = mesh%ni
        A%m = mesh%ni

        ! nombre de coefficients non nuls
        call nCoeffMat(A, mesh)

        ! allocation mémoire
        allocate(A%tmat(A%nCoeff))
        allocate(A%jposi(mesh%ni + 1))
        allocate(A%jvcell(A%nCoeff))

        ! construction jposi
        call jposiMat(A, mesh)

        ! construction jvcell
        call jvcellMat(A, mesh)

        ! construction tmat
        call tmatMat(A, mesh, d1, d2, d3)
    end subroutine



    ! =======================================================================================================
    ! FONCTIONS RELATIVES À L'ASSENBLAGE DU SCOND MEMBRE l
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! l(phi) = F - B     ici on construit F
    ! -------------------------------------------------------------------------------------------------------
    subroutine buildF(mesh, fc, F)
        implicit none

        ! paramètres
        type(Maillage), intent(in) :: mesh
        real(rp), external :: fc
        real(rp), dimension(:), intent(out) :: F

        ! variables locales
        integer :: k, i
        integer, dimension(3) :: s, slocal
        real(rp) :: x1, y1, x2, y2, x3, y3, p, aire, a, b, c

        F = 0.0_rp

        do k = 1, mesh%ne
            s(1) = mesh%connect(k, 1)
            s(2) = mesh%connect(k, 2)
            s(3) = mesh%connect(k, 3)

            ! coordonnées du premier point
            x1 = mesh%coord(mesh%connect(k, 1), 1)
            y1 = mesh%coord(mesh%connect(k, 1), 2)

            ! coordonnées du deuxième point
            x2 = mesh%coord(mesh%connect(k, 2), 1)
            y2 = mesh%coord(mesh%connect(k, 2), 2)

            ! coordonnées du troisième point
            x3 = mesh%coord(mesh%connect(k, 3), 1)
            y3 = mesh%coord(mesh%connect(k, 3), 2)

            ! aire du triangle formule de héron
            call norme2((/x1, y1/) - (/x2, y2/), a)
            call norme2((/x2, y2/) - (/x3, y3/), b)
            call norme2((/x3, y3/) - (/x1, y1/), c)
            p = (a + b + c) / 2
            aire = sqrt(p * (p - a) * (p - b) * (p - c))

            do i = 1, 3
                if (mesh%bord(s(i)) < 0) then
                    slocal(i) = -mesh%bord(s(i))
                    F(slocal(i)) = F(slocal(i)) + aire * fc(mesh%coord(s(i), 1), mesh%coord(s(i), 2)) / 3
                end if
            end do
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! l(phi) = F - B     ici on construit B
    ! -------------------------------------------------------------------------------------------------------
    subroutine buildB(mesh, t0, B)
        implicit none

        ! paramètres
        type(Maillage), intent(in) :: mesh
        real(rp), external :: t0
        real(rp), dimension(:), intent(out) :: B

        ! variables locales
        real(rp), dimension(mesh%nb, mesh%ni) :: A
        real(rp), dimension(3, 3) :: Aloc
        integer :: k, i, j, ilocal, jlocal
        integer, dimension(3) :: s

        real(rp), dimension(2) :: d1 = (/-1.0_rp, -1.0_rp/)
        real(rp), dimension(2) :: d2 = (/ 1.0_rp,  0.0_rp/)
        real(rp), dimension(2) :: d3 = (/ 0.0_rp,  1.0_rp/)

        real(rp), dimension(mesh%nb) :: cb

        !do k = 1, mesh%ne
        !    s(1) = mesh%connect(k, 1)
        !    s(2) = mesh%connect(k, 2)
        !    s(3) = mesh%connect(k, 3)

        !    Aloc = AmatLoc(k, mesh, d1, d2, d3)
        !    ilocal = 1
        !    do i = 1, 3
        !        ! colonne j
        !        do j = 1, 3
        !            if ((mesh%bord(s(i)) > 0) .and. (mesh%bord(s(j)) < 0)) then
        !                jlocal = -mesh%bord(s(j))
        !                !call ajout(ilocal, jlocal, Aloc(i, j), A)
        !                A(ilocal, jlocal) = A(ilocal, jlocal) + Aloc(i, j)
        !                ilocal = ilocal + 1
        !            end if
        !        end do
        !    end do
        !end do

        ilocal = 1
        do i = 1, mesh%ns
            if (mesh%bord(i) > 0) then
                cb(ilocal) = t0(mesh%coord(i, 1), mesh%coord(i, 2))
                ilocal = ilocal + 1
            end if
        end do

        B = 0.0_rp
        A = 0.0_rp



        do k = 1, mesh%ne
            s(1) = mesh%connect(k, 1)
            s(2) = mesh%connect(k, 2)
            s(3) = mesh%connect(k, 3)

            Aloc = AmatLoc(k, mesh, d1, d2, d3)

            if (mesh%bord(s(1)) >= 1) then
                if (mesh%bord(s(2)) < 0) then
                    A(mesh%bord(s(1)), -mesh%bord(s(2))) = A(mesh%bord(s(1)), -mesh%bord(s(2))) + Aloc(1, 2)&
                        * cb(mesh%bord(s(1)))
                end if
                if (mesh%bord(s(3)) < 0) then
                    A(mesh%bord(s(1)), -mesh%bord(s(3))) = A(mesh%bord(s(1)), -mesh%bord(s(3))) + Aloc(1, 3)&
                        * cb(mesh%bord(s(1)))
                end if
            end if

            if (mesh%bord(s(2)) >= 1) then
                if (mesh%bord(s(1)) < 0) then
                    A(mesh%bord(s(2)), -mesh%bord(s(1))) = A(mesh%bord(s(2)), -mesh%bord(s(1))) + Aloc(2, 1)&
                        * cb(mesh%bord(s(2)))
                end if
                if (mesh%bord(s(3)) < 0) then
                    A(mesh%bord(s(2)), -mesh%bord(s(3))) = A(mesh%bord(s(2)), -mesh%bord(s(3))) + Aloc(2, 3)&
                        * cb(mesh%bord(s(2)))
                end if
            end if

            if (mesh%bord(s(3)) >= 1) then
                if (mesh%bord(s(2)) < 0) then
                    A(mesh%bord(s(3)), -mesh%bord(s(2))) = A(mesh%bord(s(3)), -mesh%bord(s(2))) + Aloc(3, 2)&
                        * cb(mesh%bord(s(3)))
                end if
                if (mesh%bord(s(1)) < 0) then
                    A(mesh%bord(s(3)), -mesh%bord(s(1))) = A(mesh%bord(s(3)), -mesh%bord(s(1))) + Aloc(3, 1)&
                        * cb(mesh%bord(s(3)))
                end if
            end if
        end do

        do j = 1, mesh%ni
            B(j) = sum(A(:, j))
        end do
    end subroutine





    ! -------------------------------------------------------------------------------------------------------
    ! construction du second membre l
    ! -------------------------------------------------------------------------------------------------------
    subroutine buildL(mesh, fc, t0, L)
        implicit none

        ! paramètres
        type(Maillage), intent(in) :: mesh
        real(rp), external :: fc
        real(rp), external :: t0
        real(rp), dimension(:), allocatable :: L

        ! variables locales
        real(rp), dimension(:), allocatable :: F, B

        allocate(F(mesh%ni), B(mesh%ni), L(mesh%ni))
        call buildF(mesh, fc, F)
        call buildB(mesh, t0, B)

        L = F - B

        deallocate(F, B)
    end subroutine



    ! =======================================================================================================
    ! MANIPULATION SOLUTION
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! ajouter les termes de bord
    ! -------------------------------------------------------------------------------------------------------
    subroutine edge(mesh, t0, u, sol)
        implicit none

        ! paramètres
        type(Maillage), intent(in) :: mesh
        real(rp), external :: t0
        real(rp), dimension(:), intent(in) :: u
        real(rp), dimension(:), allocatable, intent(out) :: sol

        ! variables locales
        integer :: is, temp, isLocal

        allocate(sol(mesh%ns))

        do is = 1, mesh%ns
            temp = mesh%bord(is)
            if (temp > 0) then
                sol(is) = t0(mesh%coord(is, 1), mesh%coord(is, 2))
            else
                isLocal = -temp
                sol(is) = u(isLocal)
            end if
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! construction solution exacte
    ! -------------------------------------------------------------------------------------------------------
    subroutine buildSol(mesh, uex, sol)
        implicit none

        ! paramètres
        type(Maillage), intent(in) :: mesh
        real(rp), external :: uex
        real(rp), dimension(:), allocatable, intent(out) :: sol

        ! variables locales
        integer :: is, temp, isLocal

        allocate(sol(mesh%ns))

        do is = 1, mesh%ns
            temp = mesh%bord(is)
            isLocal = -temp
            sol(is) = uex(mesh%coord(is, 1), mesh%coord(is, 2))
        end do
    end subroutine

END MODULE elliptique
