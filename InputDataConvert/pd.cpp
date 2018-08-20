// this file is distributed under
// GPL license
#include <iostream>
#include <math_h/vectors.h>
#include <gnuplot_wrap.h>
#include <Genetic/fit.h>
#include <Genetic/uncertainties.h>
#include <Genetic/initialconditions.h>
#include <Experiment/experiment_conv.h>
#include <Kinematics/particles.h>
using namespace std;
using namespace MathTemplates;
using namespace Genetic;
using namespace GnuplotWrap;
int main()
{
    Plot th_vs_th_l("pd-th-th-lab"),t_vs_th("pd-t-th-cm-p"),cs_vs_t("pd-cs-t"),cs_vs_t_lin("pd-cs-t-lin");
    cs_vs_t<<"set key on"<<"set log y">>"unset log y";
    cs_vs_t_lin<<"set key on";
    const vector<string> files{"pb1271","pb1455","pb1459","pb1463","pb1569","pb1686"};
    FitPoints1D pointstofit;
    for(const auto&f:files){
        const auto data=Plotter::Instance().GetPoints<double,value<>>("pd/"+f);
        cs_vs_t.Points(data,f);
        cs_vs_t_lin.Points(data,f);
        for(const auto&p:data)pointstofit.push_back(p);
    }
    Fit<DifferentialMutations<>,ChiSquare,UncertaintiesEstimation> fit(pointstofit,[](const ParamSet&X,const ParamSet&P){
        return exp(P[0]+P[1]*X[0]+P[2]*X[0]*X[0]+P[3]*X[0]*X[0]*X[0]);
    });
    fit.SetUncertaintyCalcDeltas({0.001,0.001,0.001});
    fit.Init(500,make_shared<InitialDistributions>()
            <<make_shared<DistribGauss>(0.,10.)
            <<make_shared<DistribGauss>(-30.,20.)
            <<make_shared<DistribGauss>(0.,1.)
            <<make_shared<DistribGauss>(0.,1.)
    );
    while (!fit.AbsoluteOptimalityExitCondition(0.00001)) {
        fit.Iterate();
        cout << fit.iteration_count() << " iterations; "
             << fit.Optimality() << " < Chi^2 < "
             << fit.Optimality(fit.PopulationSize() - 1)
             << "           \r";
    }
    const auto chisq=fit.Optimality() / (fit.Points().size() - fit.ParamCount());
    cout << endl;
    cout << "Chi^2 divided by degrees of freedom = " << chisq << endl;
    cout << "Fit parameters with uncertainties" << endl;
    fit.SetUncertaintyCalcDeltas(parEq(fit.ParamCount(), 0.01));
    for (const auto &P : fit.ParametersWithUncertainties())
        cout << P << endl;
    cs_vs_t.Line(SortedPoints<>([&fit](double x){return fit({x});},ChainWithStep(0.0,0.001,0.4)),"chi^2/d="+to_string(chisq));
    cs_vs_t_lin.Line(SortedPoints<>([&fit](double x){return fit({x});},ChainWithStep(0.0,0.001,0.4)),"chi^2/d="+to_string(chisq));
    cs_vs_t<<"set xlabel 't,GeV/c'"<<"set ylabel 'sigma, ub/(GeV/c)'";
    t_vs_th<<"set xlabel 'theta_{p,CM},deg'"<<"set ylabel 't,GeV/c'"<<"set key on";
    th_vs_th_l<<"set xlabel 'theta_{p,lab},deg'"<<"set ylabel 'theta_{d,lab},deg'"<<"set key on"<<"set yrange [0:180]";

    for(double pb=p_beam_low;pb<=p_beam_hi;pb+=0.05){
        SortedPoints<> out;
        Points<> outpd;
        for(double theta_cm=0.001;theta_cm<PI();theta_cm+=0.001){
            const auto p0=lorentz_byPM(x()*pb,Particle::p().mass());
            const auto d0=lorentz_byPM(zero(),Particle::d().mass());
            const auto total=p0+d0;
            const auto beta=total.Beta();
            const auto p0_cm=p0.Transform(beta);
            const auto d0_cm=d0.Transform(beta);
            const auto final_cm=binaryDecay(total.M(),Particle::p().mass(),Particle::d().mass(),direction(theta_cm));
            const auto t=(final_cm.first.P()-p0.P()).M();
            const auto p1=final_cm.first.Transform(-beta);
            const auto d1=final_cm.second.Transform(-beta);
            if(t<0.7){
                out<<make_point(theta_cm*180/PI(),t);
                const auto angles=make_pair(direction(p1.P()).phi()*180./PI(),direction(d1.P()).phi()*180./PI());
                outpd.push_back(make_point(abs(angles.first),abs(angles.second)));
            }
        }
        if(out.size()>0)t_vs_th.Line(out,to_string(pb));
        if(outpd.size()>0)th_vs_th_l.Line(outpd,to_string(pb));
    }
}
