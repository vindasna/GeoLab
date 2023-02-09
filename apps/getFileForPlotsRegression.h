///////////////////////////////////////////////////////////////////////////////
//---------------------------- Libraries ------------------------------------//
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <omp.h>
#include <ncurses.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <chrono>
#include <omp.h>
#include <experimental/filesystem>

#include <map>

#include <boost/algorithm/string.hpp>



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int verbose = 0 ;

////////////////////////////////////////////////////////////////////////////////

void readTsvInfo( const char* tsvFilePath,
                  std::map< std::string, std::vector< std::string > >& tsvInfo,
                  int verbose ) ;

std::string getAgeSubject(
             const std::string subjectId,
             const std::map< std::string, std::vector< std::string > >& tsvInfo,
             int verbose ) ;


void saveOutMapDict(
      const char* savePath,
      const std::map< std::string, std::vector< std::vector< std::string > > >&
                                                                         outMap,
      const std::vector< std::string >& bundlesNames,
      int verbose ) ;
