#ifndef BUILD_GRAPH_H
#define BUILD_GRAPH_H

#include <string>

bool makeOverlapGraph(const std::string& fileIn, const std::string& fileOut, 
		  			  int minOverlap, int maxOverlap, bool filterKmer);

#endif
