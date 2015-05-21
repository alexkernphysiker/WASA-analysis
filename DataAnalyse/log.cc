#include <fstream>
#include <time.h>
#include <exception>
#include <mutex>
#include "log.h"
#include "config.h"
using namespace std;
LogLevel CurrentLogLevel=NoLog;
mutex logmutex;
string submsg[]={"","ERROR!!!","Warning","debug message"};
void SetLogLevel(LogLevel level){
	CurrentLogLevel=level;
}
void WriteToLog(LogLevel level,string msg){
	if(level==NoLog)
		throw exception();
	if(level<=CurrentLogLevel){
		lock_guard<mutex> Lock(logmutex);
		ofstream file;
		file.open("analysis.log",fstream::out|fstream::app);
		if(file.is_open()){
			time_t now=time(0);
			tm  tstruct;
			char buf[80];
			tstruct = *localtime(&now);
			strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
			file<<buf<<" "<<submsg[level]<<": "<<msg<<"\n";
			file.close();
		}
	}
}
Logger::Logger(string prefix){
	pref=prefix;
}
Logger::~Logger(){}
void Logger::LogMessage(LogLevel level, string msg){
	WriteToLog(level,pref+" -> "+msg);
}
