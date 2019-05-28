# This plot parses the output file of benchmark/verify

set terminal pdf size 30cm,18cm linewidth 1.0
set output "degree.pdf"

set grid xtics ytics

set key top left

set title 'R-MAT generator degree distribution'
set xlabel 'degree'
set ylabel 'number of nodes'

set logscale xy

plot 'degree-data.txt' with lines
