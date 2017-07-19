#!/bin/sh
#
set -e

# Converts all the TOPO precursor files to NetCDF for easy perusal
#

DIR=~/modele_input/origin


giss2nc $DIR/Z2MX2M.NGDC Z2MX2M.NGDC.nc --endian=big \
    --names=FOCEN2,ZETOP2

giss2nc $DIR/Z10MX10M Z10MX10M.nc --endian=big \
    --names=_,FLAKES

giss2nc $DIR/ZICEHXH ZICEHXH.nc --endian=big \
    --names=dZGICH,FGICEH,ZSOLDH

giss2nc $DIR/ZNGDC1 ZNGDC1.nc --endian=big \
    --names=_,_,_,FCONT1,FGICE1
