CC = g++ -std=c++11
NVCC = nvcc
CFLAGS  = -Wall -O3 -std=c++11
NVCCFLAGS = -O3 -std=c++11 -w
LDFLAGS = -lpthread

all: CONSENT

CONSENT: alignmentPiles.o reverseComplement.o kMersProcessing.o CONSENT.o DBG.o main.o
	$(NVCC) -o bin/CONSENT src/main.o src/kMersProcessing.o src/reverseComplement.o src/alignmentPiles.o src/CONSENT.o src/DBG.o BMAN-GPU/bmean.o BMAN-GPU/utils.o BMAN-GPU/BOA/align_lpo2.o  BMAN-GPU/BOA/align_lpo_po2.o  BMAN-GPU/BOA/align_score.o  BMAN-GPU/BOA/black_flag.o  BMAN-GPU/BOA/buildup_lpo.o  BMAN-GPU/BOA/create_seq.o  BMAN-GPU/BOA/fasta_format.o  BMAN-GPU/BOA/heaviest_bundle.o  BMAN-GPU/BOA/lpo_format.o  BMAN-GPU/BOA/lpo.o   BMAN-GPU/BOA/msa_format.o  BMAN-GPU/BOA/numeric_data.o  BMAN-GPU/BOA/remove_bundle.o  BMAN-GPU/BOA/seq_util.o  BMAN-GPU/BOA/stringptr.o BMAN-GPU/BOA_GPU/*.o BMAN-GPU/Complete-Striped-Smith-Waterman-Library/src/*.o $(LDFLAGS)

CONSENT.o: src/CONSENT.cpp src/CONSENT.h src/alignmentPiles.h src/kMersProcessing.h src/DBG.h
	$(NVCC) -o src/CONSENT.o -x cu -c src/CONSENT.cpp $(NVCCFLAGS) -IBMAN-GPU/BOA/ -IBMAN-GPU/BOA_GPU/

reverseComplement.o: src/reverseComplement.cpp
	$(CC) -o src/reverseComplement.o -c src/reverseComplement.cpp $(CFLAGS)

alignmentPiles.o: src/alignmentPiles.cpp src/Alignment.h src/reverseComplement.h
	$(CC) -o src/alignmentPiles.o -c src/alignmentPiles.cpp $(CFLAGS)

kMersProcessing.o: src/kMersProcessing.cpp
	$(CC) -o src/kMersProcessing.o -c src/kMersProcessing.cpp $(CFLAGS)


DBG.o: src/DBG.cpp src/reverseComplement.h
	$(CC) -o src/DBG.o -c src/DBG.cpp $(CFLAGS)

#BMAN-GPU.o: BMAN-GPU/bmean.cpp
#	$(CC) -o BMAN-GPU/bmean.o -c BMAN-GPU/bmean.cpp $(CFLAGS) -IBMAN-GPU/BOA/

utils.o: BMAN-GPU/utils.cpp
	$(CC) -o BMAN-GPU/utils.o -c BMAN-GPU/utils.cpp $(CFLAGS)

main.o: src/main.cpp src/CONSENT.h
	$(CC) -o src/main.o -c src/main.cpp

clean:
	rm src/*.o bin/CONSENT
