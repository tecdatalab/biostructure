#!/bin/sh
if [ $# -ne 1 ]; then
    echo 1>&2 Usage: $0 code.pdb
    exit 127
fi

#store arguments
receptor=$1
rec_name=`basename $receptor .pdb`

# convert to pdb.ms
rec_ms=$rec_name".pdb.ms"
./mark_sur $receptor $rec_ms


#get cp; run GETPOINTS
echo "Calculating surfaces ..."
./GETPOINTS_binary -pdb $rec_ms -smooth 0.35 -cut 1e-04

#get ZINV; run LZD
echo "Calculating Zernike ..."
rec_cp=$rec_name"_cp.txt"
rec_gts=$rec_name".gts"
rec_inv=$rec_name"_01.inv"
./LZD32 -g $rec_gts -c $rec_cp -o $rec_name -dim 161 -rad 6.0 -ord 9
rm *.dx *.grid vecCP.txt *.gts *_cp.txt *.ms
