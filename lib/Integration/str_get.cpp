// this file is distributed under 
// GPL v 3.0 license
#include <string>
#include <sstream>
#include <unistd.h>
#include "str_get.h"
using namespace std;
string ENV(string name){
	stringbuf buffer;
	ostream(&buffer)<<getenv(name.c_str());
	printf("%s: %s\n",name.c_str(),buffer.str().c_str());
}
void CD(string name){
	chdir(getenv(name.c_str()));
}
