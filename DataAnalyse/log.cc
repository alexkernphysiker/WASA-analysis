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
void Logger::AddSubprefix(string s){
	pref=pref+"=>"+s;
}
void Logger::LogMessage(LogLevel level, string msg){
	WriteToLog(level,pref+" -> "+msg);
}
Logger::SubLog Logger::getSubLog(string s){
	return Logger::getSubLog(this,s);
}
Logger::SubLog::SubLog(Logger* master, string s){
	m_master=master;
	sub_prefix=s;
	lvl=LogDebug;
	m_master->LogMessage(LogDebug,"->"+sub_prefix+" START");
}
Logger::SubLog::~SubLog(){
	m_master->LogMessage(LogDebug,"->"+sub_prefix+" END");
}
Logger::SubLog& Logger::SubLog::Message(LogLevel level, string msg){
	m_master->LogMessage(level,"->"+sub_prefix+" : "+msg);
	return *this;
}

Logger::SubLog& Logger::SubLog::operator<<(string msg){
	return Message(LogDebug,msg);
}
