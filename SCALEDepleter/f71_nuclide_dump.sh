#!/bin/bash

obiwan view -format=csv -units=atom -prec=$2 -idform='{:ee}-{:A}{:m}' "$1" | grep $3
