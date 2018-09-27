/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dedup.h
 * Author: ray
 *
 * Created on December 10, 2015, 1:37 PM
 */

#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <list>
#include <map>
#include <set>


#ifndef DEDUP_H
#define DEDUP_H

using namespace std;

#define MAXLINE 1024
//each side 24, total 48
#define TAGSIZE 24
#define TAGOFFSET 0

#define MAXN 0

void fnCountdup(string sLeft, string sRight, string sOut, int nTagSize, int nOffsetLeft, int nOffsetRight);
inline void * fnFopen(const char * szFileName, const char * szMode);
inline char * fnFgets(char * szbuffer, size_t nLen, const void * file);
inline bool fnCheckIsGZip(const char * szFile);

#endif /* DEDUP_H */

