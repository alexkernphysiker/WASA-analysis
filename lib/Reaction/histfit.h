// this file is distributed under 
// GPL v 3.0 license
#ifndef rynlxsne
#define rynlxsne
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <functional>
#include <fit.h>
#include <paramfunc.h>
#include <filter.h>
#include <initialconditions.h>
#include <equation.h>
#include <gethist.h>
void AnalyseMMSpectra(hist::point&out_bin,const hist&data,const std::vector<hist>&MC, Genetic::RANDOM&engine,std::function<void(hist&)>handle);
#endif 