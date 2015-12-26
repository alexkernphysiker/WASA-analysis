// this file is distributed under 
// MIT license
#ifndef DIQFUKBJVSPNWOVE
#define DIQFUKBJVSPNWOVE
#include <string>
#include <functional>
#include "../config.h"
#include "log.h"
#include "analysis.h"
void SetAnalysis(std::function<IAnalysis*()>analysis);
void InitLog(LogLevel level,char* type);
const std::string DataFiles=std::string("../")+DataDirectoryName+"/";
#endif
