#include <string>
enum LogLevel{NoLog=0,LogError=1,LogWarning=2,LogDebug=3};
class Logger{
public:
	Logger(std::string prefix);
	virtual ~Logger();
	virtual LogMessage(LogLevel level,std::string msg);
private:
	std::string pref;
};