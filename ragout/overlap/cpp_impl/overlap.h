#ifndef OVERLAP_H
#define OVERLAP_H

#include "fasta.h"
#include <vector>

struct Overlap
{
	Overlap(FastaRecord* prevContig, FastaRecord* nextContig, int size):
		prevContig(prevContig), nextContig(nextContig), size(size) {}

	FastaRecord* prevContig;
	FastaRecord* nextContig;
	int size;
};

std::vector<Overlap> getOverlaps(std::vector<FastaRecord>& contigs, 
								 int minOverlap, int maxOverlap);

#endif
