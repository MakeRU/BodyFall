#!/bin/sh
for file in `ls ../Data/`
do
    gnuplot -e  "out_file='../Img/${file}.png'; in_file='../Data/${file}'" MyPlotDots.gpl
#    gnuplot -e  "out_file='Img/Vel${file}.png'; in_file='Data_dots/${file}'" MyPlotVel.gpl
done
