#!/bin/bash
obiwan view -format=csv -units=atom -prec=8 -idform='{:ee}-{:A}{:m}' "$1"
