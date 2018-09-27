/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "dedup.h"
#include "zlib.h"
#include <fstream>
bool bGZ = false;


void fnCountdup(string sLeft, string sRight, string sOut, int nTagSize, int nOffsetLeft, int nOffsetRight) {
    
    FILE  *pOut;
    void *pFile1, *pFile2;
    unsigned long nTotalSampledReads = 0;
    
    map< string, unsigned int > oMap;
    //unsigned long nTotalPairs=0;
    
        
    bGZ= fnCheckIsGZip(sLeft.c_str());
    
    if (bGZ) {
        printf("GZ file detected\n");
        bGZ = true;
    }
    
    
    pFile1 = fnFopen(sLeft.c_str(), "r");
    pFile2 = fnFopen(sRight.c_str(), "r");
    pOut = fopen(sOut.c_str(), "w");
    
    if (!pFile1) {
        cout<<"Cannot open fastq 1\n";
        return;
    }
    
    if (!pFile2) {
        cout<<"Cannot open fastq 2\n";
        return;
    }
    
    if (!pOut) {
        cout<<"Cannot open output file for writing\n";
        return;
    }

    char szBuffer1[MAXLINE];
    char szBuffer2[MAXLINE];

    bool bReady = false; // ready to read line?
    bool bQualLine = false;
    while (fnFgets(szBuffer1, MAXLINE, pFile1)) {
	if (!fnFgets(szBuffer2, MAXLINE, pFile2)) {
            cout<<"Error:Two fastq files not in sync.\n";
            return;
        }
        
        //            printf("buffer1 len %d, buffer2 len %d\n" , strlen(szBuffer1), strlen(szBuffer2));
        
        if (bQualLine) {
            bQualLine = false;
            continue;
        }
                    
        if ( szBuffer1[0] == '@'  ||  szBuffer2[0] == '@' ) {
            
            if (szBuffer1[0] != szBuffer2[0]) {
                cout<<"Error:Two fastq files not in sync.\n";
                return;
            }
            
            //empty string buffer
            //memset(szBuffer1,0, MAXLINE );
            //memset(szBuffer2,0, MAXLINE );
            
            bReady = true;
            continue;
        }
                    
        if ( (szBuffer1[0]=='+' && strlen(szBuffer1)<=3 ) ||  (szBuffer2[0]=='+' && strlen(szBuffer2)<=3 ) ) {
            if (szBuffer1[0] != szBuffer2[0]) {
                cout<<"Error:Two fastq files not in sync.\n";
                return;
            }   
            bReady = false;
            bQualLine = true;
            continue;
        }
                     
        
        if (bReady) {
            bReady = false;
            string sSeq1(szBuffer1);
            string sSeq2(szBuffer2);
            
            //empty string buffer
            //memset(szBuffer1,0, MAXLINE );
            //memset(szBuffer2,0, MAXLINE );
            if (sSeq1.size() < nTagSize || sSeq2.size() <  nTagSize) { // ignore
                continue;
            }
            
            
            string sTag;
            sTag.append(sSeq1.substr(nOffsetLeft , nTagSize));
            sTag.append(sSeq2.substr(nOffsetRight, nTagSize));
            
            
            size_t nN = std::count(sTag.begin(), sTag.end(), 'N'); // COUNT Ns
            if (nN > MAXN) {
                continue;
            }
                        
            if (oMap.find(sTag) == oMap.end()) {
                oMap[sTag] = 1;
            } else {
                oMap[sTag]++;
            }
            
            nTotalSampledReads++;
        }
    }
    
    for (map< string, unsigned int >::iterator it=oMap.begin(); it!=oMap.end(); ++it) {
        fprintf (pOut, "%u\n",it->second);
    }
    
    printf( "%lu" , nTotalSampledReads );

}

inline void * fnFopen(const char * szFileName, const char * szMode) {
    if (bGZ) {
        return (gzopen(szFileName , szMode));
    } 
    
    return fopen(szFileName , szMode);
}

inline char * fnFgets(char * szbuffer, size_t nLen, const void * file) {
    if (bGZ) {
        return gzgets( (gzFile)file,  szbuffer , nLen );
    } 
    
    return fgets(szbuffer , nLen, (FILE*)(file));
}


inline bool fnCheckIsGZip(const char * szFile) {
    ifstream file(szFile, ios_base::binary | ios_base::in);
    char nByte1, nByte2;
    file.seekg(0);
    file.read(&nByte1, 1);
    file.seekg(1);
    file.read(&nByte2, 1);
    file.close();
    //cout<< nByte1 << " " << nByte2 << "\n";
    if ( (nByte1 == '\x1f' )&& (nByte2 == '\x8b')) {
        return true;
    }
    else {
        return false;
    }
}
