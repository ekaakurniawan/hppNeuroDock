python Molecule.py | grep "^ATOM" | cut -c31-54 > log
sed -n '6696,6768p' Results/ind.dlg| cut -c31-54 > ref
diff log ref
rm log
rm ref
