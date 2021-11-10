for line in `cat fname.txt`;do
	sed -i '/ATP/d' $line
	sed -i '/STR/d' $line
	sed -i '/DMP/d' $line
	sed -i '/IPT/d' $line
done
