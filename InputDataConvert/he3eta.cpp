// this file is distributed under
// GPL license
#include <iostream>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    const LinearInterpolation<>
    alpha=Points<>{
        {-0.001,0},{0.04,0},
        {0.112,0.52},{0.130,0.62},{0.145,0.74},{0.170,1.2}
    },
    beta=Points<>{
        {-0.001,-0.32},{0.04,-0.32},
        {0.112,-0.33},{0.130,-0.34},{0.145,-0.30},{0.170,-0.34}
    },
    gamma=Points<>{
        {-0.001,0},{0.04,0},
        {0.112,-0.10},{0.130,-0.18},{0.145,-0.21},{0.170,-0.29}
    };
    const auto chain=ChainWithStep(0.0,0.001,0.170);
    Plot("he3eta-alpha").Line(SortedPoints<>(alpha.func(),chain)).Points(alpha)
        << "set xlabel 'p_{eta,CM}, GeV/c'"<<"set xrange [0.1:]"
        << "set ylabel 'alpha, n.d.'"<<"set yrange [0:1.5]"
        << "set title 'Alpha parameter'";
    Plot("he3eta-beta").Line(SortedPoints<>(beta.func(),chain)).Points(beta)
        << "set xlabel 'p_{eta,CM}, GeV/c'"<<"set xrange [0.1:]"
        << "set ylabel 'beta, n.d.'"<<"set yrange [-0.5:0]"
        << "set title 'Beta parameter'";
    Plot("he3eta-gamma").Line(SortedPoints<>(gamma.func(),chain)).Points(gamma)
        << "set xlabel 'p_{eta,CM}, GeV/c'"<<"set xrange [0.1:]"
        << "set ylabel 'gamma, n.d.'"<<"set yrange [-0.3:0]"
        << "set title 'Gamma parameter'";
}

