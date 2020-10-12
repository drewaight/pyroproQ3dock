#!/usr/bin/env bash


read -p "File to change:" file_var
read -p "Chain to change:" inp_var
read -p "Change to chain:" out_var

sed -i "s/^\(ATOM.\{17\}\)$inp_var/\1$out_var/" "$file_var"
