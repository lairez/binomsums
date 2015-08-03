
SRC = src/binomsums.mpl
MLA = binomsums.mla

$(MLA) : $(SRC)
	rm -f $(MLA)
	maple -sq <<< 'march(create, "$(MLA)", 10);'
	maple -b $(MLA) -B -D BUILD_MLA -w 4 -sq $(SRC)
        


