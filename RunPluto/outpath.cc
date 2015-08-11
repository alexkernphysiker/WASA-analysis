// this file is distributed under 
// GPL v 3.0 license
string outpath;{
	stringbuf buffer;
	ostream os (&buffer); 
	os<<getenv("PLUTO_OUTPUT");
	outpath=buffer.str();
	printf("output path: %s\n",outpath.c_str());
}