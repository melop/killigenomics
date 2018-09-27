
function linkfiles {
	DIR=$1;
	mkdir -p $DIR
	ln -sf `realpath ../allelefreq`/$DIR/out_f_GFE_*_of_*.txt $DIR/
}

linkfiles RACWET.excludehybrid
linkfiles RACDRY.excludehybrid
linkfiles ORTDRY.excludehybrid
linkfiles ORTWET

