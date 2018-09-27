/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: ray
 *
 * Created on December 10, 2015, 1:33 PM
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include "dedup.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    if (argc == 4) {
	fnCountdup(argv[1] , argv[2], argv[3], TAGSIZE, TAGOFFSET, TAGOFFSET);
	cout<<"Running: ./countduppe "<< argv[1] << " " << argv[2] << " " << argv[3] << " " << TAGSIZE << " "<< TAGOFFSET << "\n";

	return 0;
    } else if (argc == 5) {
        int nTagSize = TAGSIZE;
        sscanf(argv[4] , "%i", &nTagSize );
        cout<<"Running: ./countduppe "<< argv[1] << " " << argv[2] << " " << argv[3] << " " << nTagSize << " " << TAGOFFSET << "\n";

       	fnCountdup(argv[1] , argv[2], argv[3], nTagSize, TAGOFFSET, TAGOFFSET);
	
	return 0; 
    } else if (argc == 7) {
        int nTagSize = 24;
        sscanf(argv[4] , "%i", &nTagSize);
        int nTagOffsetLeft = TAGOFFSET;
        sscanf(argv[5] , "%i", &nTagOffsetLeft);
        int nTagOffsetRight = TAGOFFSET;
        sscanf(argv[6] , "%i", &nTagOffsetRight);
        
        cout<<"Running: ./countduppe "<< argv[1] << " " << argv[2] << " " << argv[3] << " " << nTagSize << " " <<  nTagOffsetLeft << " " << nTagOffsetRight << "\n";

       	fnCountdup(argv[1] , argv[2], argv[3], nTagSize, nTagOffsetLeft, nTagOffsetRight);
	
	return 0; 
    }
    
    else {
        
        cout<<"Usage: ./countduppe 1.fastq 2.fastq out.txt [tagsize in bp = 24] [tag offset in read1 = 0] [tag offset in read 2 =0] \n";
    }
        
    return 0;
}

