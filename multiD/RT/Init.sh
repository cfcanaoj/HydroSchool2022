#!/bin/bash


dir=xydata
if [ ! -d ${dir} ]; then
    mkdir ${dir}
fi

dir=figures
if [ ! -d ${dir} ]; then
    mkdir ${dir}
fi

dir=movies
if [ ! -d ${dir} ]; then
    mkdir ${dir}
fi


