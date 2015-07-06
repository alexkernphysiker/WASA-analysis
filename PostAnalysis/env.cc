// this file is distributed under 
// GPL v 3.0 license
	string inputpath;{
		stringbuf buffer;
		ostream os (&buffer); 
		os<<getenv("WASA_OUTPUT_DATA");
		inputpath=buffer.str();
		printf("input path: %s\n",inputpath.c_str());
	}
	string outpath;{
		stringbuf buffer;
		ostream os (&buffer); 
		os<<getenv("POST_ANALYSIS_DATA");
		outpath=buffer.str();
		printf("output path: %s\n",outpath.c_str());
	}