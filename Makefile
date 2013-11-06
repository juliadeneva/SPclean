spclean: spclean.h
	gcc -o $@ spclean.c ANDclean.c DFTclean.c -lm
