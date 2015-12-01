// this file is distributed under 
// MIT license
#include <string>
#include <sstream>
#include <unistd.h>
#include <list>
#include "str_get.h"
using namespace std;
list<string> dir_cache;
void PUSHD(){
	dir_cache.push_front(string(getcwd(NULL,0)));
}
void POPD(){
	chdir(dir_cache.begin()->c_str());
	dir_cache.pop_front();
}

string ENV(string name){
	stringbuf buffer;
	ostream(&buffer)<<getenv(name.c_str());
	printf("%s: %s\n",name.c_str(),buffer.str().c_str());
	return buffer.str();
}
void CD(string name){
	chdir(getenv(name.c_str()));
}
