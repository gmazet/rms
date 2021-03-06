
SOFT=`basename $0`

case $# in
5)      year=$1 ;
        month=$2 ;
	day=$3;
	duration=$4;
        chans=$5 ;
	month=$(echo $month | awk '{printf ("%2.2d",$1)}');
	lastNdays=0;
	daybefore=32 ; # HACK
	BEGDATE=$(echo $year $month $day | awk '{printf ("%4.4d-%2.2d-%2.2dT00:00:00.000", $1,$2,$3)}');;
2)	lastNdays=$1;
        chans=$2 ;
	duration=$(echo $lastNdays | awk '{print ($1-1)*24}')
	BEGDATE=$(date '+%Y-%m-%d' -d '-'$lastNdays'days');
	year=$(echo $BEGDATE | awk -F"-" '{print $1}');
	month=$(echo $BEGDATE | awk -F"-" '{print $2}');
	daybefore=$(date '+%d' -d '-1day');;
*)      echo "($SOFT) usage : $SOFT <YYYY> <MM> <DD> <DURATION in hours> <GROUP [SH|BH|MURU]> or $SOFT <NB DAYS BACK> <GROUP [SH|BH|MURU]>";
        exit;;
esac

echo "BEGDATE=$BEGDATE"
echo "year,month=$year,$month"

# Répertoire de travail
DIR=/cea/data/outils/dase/noise_analysis/rms
# Où l'on stocke tous les résultats de PPSD et de RMS
DATADIR=/cea/dsku/SDS_muru/SDS_muru/rms


if [ $lastNdays -eq 0 ] ; then
	# Taille de la fenêtre glissante en secondes pour la rms
	WINLEN=1800
else
	# Taille de la fenêtre glissante en secondes pour la rms
	WINLEN=600
fi

PYTHONEXE=/usr/bin/python3

if [ $chans = "MURU" ] ; then
	# Détecteur RC de MURU
	FMIN=8.0
	FMAX=16.0
	RMS_UNIT=VEL
	strgrep="[HE][HP][ZNE]"
	grep "^RD\|" $DIR/stations_DASE.nslc.txt | grep "MURU" > $DIR/tmp/$chans.txt
else
	# Détecteur RSN d'alerte
	FMIN=2.0
	FMAX=6.0
	RMS_UNIT=DISP
        if [ $chans = "SH" ] ; then
                strgrep="\|[ES]H[ZNE]\||SFTF"
        else
                strgrep="[BH]H[ZNE]"
        fi
	grep "^RD\|" $DIR/stations_DASE.nslc.txt | grep -v "MURU" | egrep $strgrep  > $DIR/tmp/$chans.txt
fi
cat $DIR/tmp/$chans.txt


# On ne prend que les voies concernées
nbchan=$(wc -l $DIR/tmp/$chans.txt | awk '{print $1}')
i=1
while ( [ $i -le $nbchan ] ) ; do
#for sta in ETNF ETSF ; do
	net=$(tail -n +$i $DIR/tmp/$chans.txt | head -1 | awk -F"|" '{print $1}')
	sta=$(tail -n +$i $DIR/tmp/$chans.txt | head -1 | awk -F"|" '{print $2}')
	loccode=$(tail -n +$i $DIR/tmp/$chans.txt | head -1 | awk -F"|" '{print $3}')
	channel=$(tail -n +$i $DIR/tmp/$chans.txt | head -1 | awk -F"|" '{print $4}')
	##if [ $sta != "MEZF" ] ; then i=$(expr $i + 1);continue; fi
	echo
	echo
	echo "#####################################################################"
        echo "Process " $net $sta "$loccode" $channel

	# fréquence de Nyquist pour chaque voie
	chancode=$(echo $channel | cut -c 1-1)
	if [ $chancode = "B" -o $chancode = "S" ] ; then
		freqmax=25
	else
		freqmax=50
	fi

        echo $BEGDATE $duration $net $sta "$loccode" $channel $freqmax
	#i=$(expr $i + 1)
	#continue

	echo "----------------------------"
	echo "--- Plot spectrum"
	echo "----------------------------"
 	sleep 0.5;$PYTHONEXE $DIR/spec.py $BEGDATE $duration $net $sta "$loccode" $channel $freqmax
	
	echo "----------------------------"
	echo "--- Plot spectro+rms"
	echo "----------------------------"
	FREQS=$(echo $FMIN $FMAX | awk '{printf ("%.1f_%.1fHz", $1,$2)}')
	echo "FREQS=",$FREQS
 	sleep 0.5;$PYTHONEXE $DIR/combo2.py $BEGDATE $duration $net $sta "$loccode" $channel $FMIN $FMAX "$FREQS" "$FREQS" $freqmax $RMS_UNIT

	if [ $lastNdays -eq 0 ] ; then
		mkdir -p $DIR/figures/$year/$month
		mv $DIR/figures/$net.$sta.$loccode.*.png $DIR/figures/$year/$month/.
	else
		mkdir -p $DIR/figures/last.$lastNdays.days
		mv $DIR/figures/$net.$sta.$loccode.*.png $DIR/figures/last.$lastNdays.days/.
	fi

#echo "Presse a key..."
#read A
#exit
	date
	i=$(expr $i + 1)

done

if [ $lastNdays -eq 0 ] ; then
	$DIR/figures/build_html.sh $year $month $chans
else
	$DIR/figures/build_html.sh $lastNdays $chans
fi

