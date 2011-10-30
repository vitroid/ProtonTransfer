60x1234.ar3a: ../visualize
	../visualize --site --compose=60 1234 > $@
60x1234.ngph: ../visualize
	../visualize --graph --compose=60 1234 > $@
60x16x1234.log: 60x16x1234.sh
	sh $< > $@
#probability of reaching by the shortest paths
60x16x1234.pa: 60x16x1234.log
	python ../pathanalysis.py < $< > $@
#raw number of possible paths
60x16x1234.paraw: 60x16x1234.log
	python ../pathanalysis.py -r < $< > $@
#probability of reaching within n steps
60x16x1234.pasum: 60x16x1234.log
	python ../pathanalysis.py -s < $< > $@
