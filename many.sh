
BEGTIME="2021-01-12T00:00:00"
WINLEN=600
NBHOURS=0

NET="AM"
LOCCODE="00"
CHAN="SHZ"

NET="FR"
LOCCODE="00"
CHAN="HHZ"

OUTPUT="DISP"

FMIN=2.0
FMAX=10.0

N=23
M=200
for sta in $(cat ./resif.sta | awk -F"|" '{print $2}' | tail -n +$N | head -$M) ; do
#for sta in R9F1B ; do
	echo "##############################################################################"
	echo $sta
	echo "##############################################################################"

	./deconvolution.py $BEGTIME $NBHOURS $NET $sta $LOCCODE $CHAN $OUTPUT
	sleep 1
	./rms_in_csv.py $BEGTIME $NBHOURS $NET $sta $LOCCODE $CHAN $FMIN $FMAX $WINLEN
	sleep 1

	./ppsd_and_npz.py $BEGTIME $NBHOURS $NET $sta $LOCCODE $CHAN
	sleep 1
	rm ./rawdata/2021/*/*/$net.$sta.*.mseed
	rm ./DISP/2021/*/*/$net.$sta.*.mseed
	./combo.py $BEGTIME $NBHOURS $NET $sta $LOCCODE $CHAN $FMIN $FMAX "2.0_10.0Hz" "2.0_10.0Hz"
	sleep 1
	echo
done
