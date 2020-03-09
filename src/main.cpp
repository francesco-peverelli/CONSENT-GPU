#include "CONSENT.h"
#include<chrono>
#include<fstream>
#include<ctime>

#define GPU 1
#define NOW std::chrono::high_resolution_clock::now()

int main(int argc, char* argv[]) {

	if (argc < 2) {
		fprintf(stderr, "Usage: %s [-a alignmentFile.paf] [-s minSupportForGoodRegions] [-l minLengthForGoodRegions] [-j threadsNb] [-B maxBatch]\n\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	std::string PAFIndex, alignmentFile, readsFile, proofFile, path;
	PAFIndex = "";
	alignmentFile = "";
	readsFile = "";
	unsigned minSupport, maxSupport, maxMSA, windowSize, nbThreads, opt, merSize, commonKMers, minAnchors, solidThresh, windowOverlap, nbReads, maxBatch;

	minSupport = 4;
	maxSupport = 1000;
	maxMSA = 150;
	windowSize = 500;
	merSize = 9;
	commonKMers = 8;
	minAnchors = 10;
	solidThresh = 4;
	windowOverlap = 50;
	nbThreads = 1;
	nbReads = 0;
	maxBatch = 5000000;


	while ((opt = getopt(argc, argv, "a:B:A:d:k:s:S:M:l:f:e:p:c:m:j:w:m:r:R:n:i:")) != -1) {
        switch (opt) {
        	case 'i':
        		PAFIndex = optarg;
        		break;
			case 'a':
				alignmentFile = optarg;
				break;
			case 'B':
				maxBatch = atoi(optarg);
				break;
			case 's':
				minSupport = atoi(optarg);
				break;
			case 'S':
				maxSupport = atoi(optarg);
				break;
			case 'M':
				maxMSA = atoi(optarg);
				break;
			case 'l':
				windowSize = atoi(optarg);
				break;
			case 'k':
				merSize = atoi(optarg);
				break;
			case 'c':
				commonKMers = atoi(optarg);
				break;
			case 'A':
				minAnchors = atoi(optarg);
				break;
			case 'f':
				solidThresh = atoi(optarg);
				break;
			case 'm':
				windowOverlap = atoi(optarg);
				break;
			case 'r':
				readsFile = optarg;
				break;
			case 'R':
				proofFile = optarg;
				break;
			case 'p':
				path = optarg;
				path +=  + "/BMEAN/BOA/blosum80.mat";
				break;
			case 'n':
				nbReads = atoi(optarg);
				break;
			case 'j':
				nbThreads = atoi(optarg);
				break;
			default: /* '?' */
				fprintf(stderr, "Usage: %s [-a alignmentFile.paf] [-k merSize] [-s minSupportForGoodRegions] [-l minLengthForGoodRegions] [-f freqThresholdForKMers] [-e maxError] [-p freqThresholdForKPersFreqs] [-c freqThresholdForKPersCons] [-m mode (0 for regions, 1 for cluster)] [-j threadsNb] [-B maxBatch]\n\n", argv[0]);
				exit(EXIT_FAILURE);
        }
    }
	auto corr_s = NOW;

#if GPU 
	
	runCorrection_gpu(PAFIndex, alignmentFile, minSupport, maxSupport, windowSize, merSize, commonKMers, minAnchors, solidThresh, windowOverlap, nbThreads, readsFile, proofFile, maxMSA, maxBatch, path);

#endif
#if (GPU == 0)

	runCorrection(PAFIndex, alignmentFile, minSupport, maxSupport, windowSize, merSize, commonKMers, minAnchors, solidThresh, windowOverlap, nbThreads, readsFile, proofFile, maxMSA, path);

#endif
	auto corr_e = NOW;

	auto tot_time = std::chrono::duration_cast<std::chrono::milliseconds>(corr_e - corr_s);

	std::ofstream outfile;
	outfile.open("c_timings.txt", std::ios_base::app);
		
	std::time_t st_time = std::chrono::system_clock::to_time_t(corr_s);
	outfile << "Timed " << std::ctime(&st_time) << ": " << tot_time.count() << " milliseconds\n";	

	return EXIT_SUCCESS;
}
