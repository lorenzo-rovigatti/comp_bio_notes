#!/bin/bash

LC_ALL=C

for f in $(ls *agr)
do
	base=${f%.*}
	if [ $base.png -ot $f ]
	then
		LANG= xmgrace -hardcopy -hdevice EPS $f -printfile $base.eps
		epstopdf $base.eps --outfile=$base.pdf
		convert $base.pdf $base.png
		rm $base.eps $base.pdf
	fi
done
