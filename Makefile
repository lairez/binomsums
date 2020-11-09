SRC = $(wildcard src/*.mpl)
MLA = binomsums.mla

all: $(MLA)

$(MLA) : $(SRC)
	rm -f $(MLA)
	echo 'march(create, "$(MLA)", 10);' | maple -sq
	maple -b $(MLA) -B -I src -D BUILD_MLA -w 4 -sq src/pack.mpl 

test : all force
	mkdir -p log
	for f in test/*.test.mpl ;      do \
		echo Testing $$f  ; \
		maple -b $(MLA) -B -i test/init.mpl -s $$f | tee log/$$(basename $$f .tst).log | grep FAILED ;\
		true; \
	done

force :
	true
