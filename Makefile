# -----------------------------------------------------------------------------------------------------------
#  variables
# -----------------------------------------------------------------------------------------------------------
PROJET  = elt_finis
SOURCES = main.f95 elliptique.f95 csr.f95 math.f95 sdd.f95 donnees.f95
MODULES =          elliptique.mod csr.mod math.mod sdd.mod donnees.mod
OBJECTS = math.o   elliptique.o   csr.o   main.o   sdd.o   donnees.o
CC      = gfortran -fbounds-check
EXEC    = truc
OTHER   =
PDF     = autre/cpp_laplacien1D.pdf



# -----------------------------------------------------------------------------------------------------------
#  compilation rules. Pour une liste des dépendences : gfortran -M -cpp fichier.f95 (après avoir créer .o)
# -----------------------------------------------------------------------------------------------------------
# exécutable
$(EXEC) : $(OBJECTS) maillage/carre.poly
	make -C maillage > /dev/null
	$(CC) $(OBJECTS) -o $(EXEC)

# programme principal
main.o : main.f95 $(MODULES)
	$(CC) -c main.f95 -o main.o

# modules
math.o math.mod : math.f95
	$(CC) -c math.f95 -o math.o

csr.o csr.mod : csr.f95 math.mod
	$(CC) -c csr.f95 -o csr.o

elliptique.o elliptique.mod : elliptique.f95 math.mod sdd.mod csr.mod
	$(CC) -c elliptique.f95 -o elliptique.o

donnees.o donnees.mod : donnees.f95 math.mod
	$(CC) -c donnees.f95 -o donnees.o

sdd.o sdd.mod : sdd.f95 math.mod
	$(CC) -c sdd.f95 -o sdd.o



# -----------------------------------------------------------------------------------------------------------
# phony targets
# -----------------------------------------------------------------------------------------------------------
# exécuter l'exécutable (utile pour utiliser F5 dans vim)
.PHONY : make
make :
	@./$(EXEC)

# supprimer les fichiers objet et l'exécutable s'ils existent
.PHONY : clean
clean :
	rm -f $(EXEC) $(OBJECTS) $(MODULES)

# effacer le contenu des dossiers d'entrees et de sorties
.PHONY : del
del :
	rm -f sorties/* graphes/*

# ouvrir les fichiers du projet dans des onglets de vim
.PHONY : open
open :
	@vim -p $(SOURCES) Makefile Plot.gnu $(OTHER)

# tout compiler et lancer gdb (segmentation fault)
.PHONY : gdb
gdb :
	$(CC) -g $(SOURCES) -o $(EXEC) -lm && gdb ./$(EXEC)

# clean et tarer le dossier
.PHONY : tar
tar :
	make clean
	rm -f .*.swp
	tar -zcvf ../$(PROJET).tar.gz ../$(PROJET)

# sauvegarder ancienne version
.PHONY : save
save :
	make clean
	rm -rf ../old_$(PROJET)
	cp -r ../$(PROJET) ../old_$(PROJET)

.PHONY : pdf
pdf :
	xdg-open $(PDF)

#
.PHONY : clean
coffe :
	@echo "  (\n   )\n c[]"

.PHONY : show
show :
	make -C maillage show
