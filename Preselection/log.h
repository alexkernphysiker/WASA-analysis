// this file is distributed under 
// MIT license
#ifndef YGCXMVEIRHRTKMLS
#define YGCXMVEIRHRTKMLS
#include <string>
enum LogLevel{NoLog=0,LogError=1,LogWarning=2,LogDebug=3};
class Logger{
public:
	class SubLog{
	public:
		SubLog(Logger*master,LogLevel l);
		virtual ~SubLog();
		SubLog&operator<<(std::string msg);
	private:
		LogLevel lvl;
		Logger *m_master;
	};
	friend class SubLog;
	Logger();
	virtual ~Logger();
	void AddLogSubprefix(std::string s);
	SubLog Log(LogLevel l=LogWarning);
protected:
	void LogMessage(LogLevel level,std::string msg);
private:
	std::string pref;
};
#endif