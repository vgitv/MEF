MODULE sdd

    use math

    implicit none

    type :: Cell
        integer :: elt
        type(Cell), pointer :: suiv
    end type Cell

    type :: Liste
        type(Cell), pointer :: premier
        integer :: taille
    end type Liste

contains

    ! -------------------------------------------------------------------------------------------------------
    ! constructeur
    ! -------------------------------------------------------------------------------------------------------
    subroutine newListe(l)
        implicit none

        ! paramètres
        type(Liste) :: l

        l%premier => NULL()
        l%taille = 0
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! destructeur
    ! -------------------------------------------------------------------------------------------------------
    subroutine rmListe(l)
        implicit none

        ! paramètres
        type(Liste) :: l

        ! variables locales
        type(Cell), pointer :: p, pmem

        p => l%premier

        do while (associated(p))
            pmem => p%suiv
            deallocate(p)
            p => pmem
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! insérer un élément après une cellule pointée par p dans la liste l
    ! -------------------------------------------------------------------------------------------------------
    subroutine insertListe(l, p, elt)
        implicit none

        ! paramètres
        type(Liste), intent(inout) :: l
        type(Cell), pointer :: p
        integer, intent(in) :: elt

        ! variables locales
        type(Cell), pointer :: cellule

        allocate(cellule)

        cellule%elt = elt

        if (associated(p)) then
            if (associated(p%suiv)) then
                ! insertion en milieu de liste
                cellule%suiv => p%suiv
                p%suiv => cellule
            else
                ! insertion en fin de liste
                cellule%suiv => NULL()
                p%suiv => cellule
            end if
        else
            ! insertion en tête
            cellule%suiv => l%premier
            l%premier => cellule
        end if

        l%taille = l%taille + 1
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! affiche la liste
    ! -------------------------------------------------------------------------------------------------------
    subroutine printListe(l)
        implicit none

        ! paramètres
        type(Liste), intent(in) :: l

        ! variables locales
        type(Cell), pointer :: p

        p => l%premier

        do while (associated(p))
            print *, p%elt
            p => p%suiv
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! insère elt dans une liste triée à la bonne place (si non présent) et renvoie true si l'elt à été ajouté
    ! -------------------------------------------------------------------------------------------------------
    function insertTrieListe(l, elt)
        implicit none

        ! paramètres
        type(Liste), intent(inout) :: l
        integer, intent(in) :: elt

        ! return
        logical :: insertTrieListe

        ! variables locales
        type(Cell), pointer :: p
        integer :: courant

        insertTrieListe = .false.
        p => l%premier
        courant = p%elt

        if (courant < elt) then
            do while (associated(p))
                if (associated(p%suiv)) then
                    if (elt == p%suiv%elt) then
                        exit
                    else if (elt < p%suiv%elt) then
                        call insertListe(l, p, elt)
                        insertTrieListe = .true.
                        exit
                    end if
                else if (p%elt < elt) then
                    call insertListe(l, p, elt)
                    insertTrieListe = .true.
                    exit
                end if
                p => p%suiv
            end do
        else
            call insertListe(l, NULL(), elt)
            insertTrieListe = .true.
        end if
    end function

END MODULE sdd
