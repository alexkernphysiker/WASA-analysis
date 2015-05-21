#ifndef YGCXMVEIRHRTKMLS
#define YGCXMVEIRHRTKMLS
#include <string>
enum LogLevel{NoLog=0,LogError=1,LogWarning=2,LogDebug=3};
class Logger{
public:
	class SubLog{
		SubLog(Logger*master,std::string s);
		virtual ~SubLog();
		SubLog&operator<<(std::string msg);
		SubLog&Message(LogLevel level,std::string msg);
	private:
		LogLevel lvl;
		std::string sub_prefix;
		Logger *m_master;
	};
	friend class SubLog;
	Logger(std::string prefix);
	virtual ~Logger();
	void AddSubprefix(std::string s);
	SubLog getSubLog(std::string s);
protected:
	void LogMessage(LogLevel level,std::string msg);
private:
	std::string pref;
};
#endif