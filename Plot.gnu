# réinitialisation des paramètres
reset

# fichier de sortie
set term postscript eps size 3.5,2.62 enhanced color font 'Helvetica,12'
set output "graphes/u.eps"
set encoding utf8

# paramètres
set title "Équation de Poisson 2D"
set grid
set xlabel "x"
set ylabel "y"
#set size ratio -1 # repère orthonormé
#set xrange[-2:2]
#set yrange[-1.5:1.5]


# tracé
splot "sorties/sol.dat"        u 1:2:3  with lines  lc rgb "#FF0000" lw 1 title "Solution approchée",\
      "sorties/uOmega.dat"     u 1:2:3  with lines  lc rgb "#008000" lw 1 title "Solution exacte"

# affichage écran
set term wxt enhanced
replot
