// this file is distributed under
// GPL license
#include <fstream>
#include "he3eta.h"
using namespace MathTemplates;
using namespace std;
const Reaction &he3eta()
{
    static Reaction main_react(Particle::p(), Particle::d(), {Particle::he3(), Particle::eta()});
    return main_react;
}
const LinearInterpolation<value<double>> &he3eta_sigma()
{
    static LinearInterpolation<value<double>> res;
    if (res.size() == 0) {
        ifstream file("he3eta.txt");
        while (file) {
            double x, y, dy;
            file >> x >> y >> dy;
            res << point<value<double>>(x, {y, dy});
        }
        file.close();
    }
    return res;
}
