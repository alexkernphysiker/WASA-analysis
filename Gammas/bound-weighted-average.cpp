// this file is distributed under
// GPL license
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <math_h/sigma3.h>
#include <Genetic/fit.h>
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include <Parameters/parameters.h>
#include <Parameters/systematic.h>
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "bound-weighted-average");
    const ext_hist<2> data1 = Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("curve_3he_2gamma");
    const value<> branching_ratio1{0.393,0.003};
    const ext_hist<2> data2 = Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("curve_3he_6gamma");
    const value<> branching_ratio2{0.322,0.003};

    const hist<> avr=hist_avr(wrap_hist(data1)/branching_ratio1,wrap_hist(data2)/branching_ratio2);
    Plot("Average",5)
        .Hist(avr)
            <<"set key left top">>"set key right top"
            << "set xlabel 'Q, MeV'" << "set key on"<<"set xrange [-70:10]"
            << "set ylabel 'Normalized events, nb'" << "set yrange [0:70]"
            << "set title 'Weighted average'"<<"set key right top";
}
