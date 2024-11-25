/**
 * \file build_shapes.cpp
 * \brief This file contains functions to create panels of different
 * shapes, and also to compute their normal vectors (in terms of angles).
 *
 * This file is a part of the simpleTBEM library.
 */
#include "math.h"
#include "parametrized_line.hpp"
#include "parametrized_circular_arc.hpp"

// define shorthand for time benchmarking tools, complex data type and immaginary unit
typedef std::complex<double> complex_t;
double eps = 1E-10;
double pififths = .2*M_PI;
double piquarters = .25*M_PI;
double twopififths = .4*M_PI;
double pihalves = .5*M_PI;
double threepiquarters = .75*M_PI;

PanelVector buildCircle(double radius, unsigned &numpanels) {
    ParametrizedCircularArc curve(Eigen::Vector2d(0,0),radius,0,2*M_PI);
    return curve.split(numpanels);
}

double normalCircle(double x1, double x2) {
    return atan2(x2, x1);
}

PanelVector buildSquare(double halfsidelength, unsigned &numpanels) {
    int numpanelseach = (numpanels + 3) / 4;
    numpanels = 4 * numpanelseach; // round up

    // Corner points
    Eigen::RowVectorXd x1(2); x1 << halfsidelength, -halfsidelength;
    Eigen::RowVectorXd x2(2); x2 << halfsidelength, halfsidelength;
    Eigen::RowVectorXd x3(2); x3 << - halfsidelength, halfsidelength;
    Eigen::RowVectorXd x4(2); x4 << - halfsidelength, - halfsidelength;

    // parametrized line segments forming the edges of the polygon
    ParametrizedLine line1(x1, x2); PanelVector line1panels = line1.split(numpanelseach);
    ParametrizedLine line2(x2, x3); PanelVector line2panels = line2.split(numpanelseach);
    ParametrizedLine line3(x3, x4); PanelVector line3panels = line3.split(numpanelseach);
    ParametrizedLine line4(x4, x1); PanelVector line4panels = line4.split(numpanelseach);

    // storing all the panels in order so that they form a polygon
    PanelVector panels;
    panels.insert(panels.end(), line1panels.begin(), line1panels.end());
    panels.insert(panels.end(), line2panels.begin(), line2panels.end());
    panels.insert(panels.end(), line3panels.begin(), line3panels.end());
    panels.insert(panels.end(), line4panels.begin(), line4panels.end());
    return panels;
}

double normalSquare(double x1, double x2) {
    double angle = atan2(x2, x1);
    int wedge_index = floor((angle + piquarters) / pihalves);
    double wedge_center = wedge_index * pihalves;
    double wedge_left = wedge_center - piquarters;
    double wedge_right = wedge_center + piquarters;

    if (angle - wedge_left < eps) return wedge_left;
    if (wedge_right - angle < eps) return wedge_right;
    return wedge_center;
}

PanelVector buildStar(double innerradius, double outerradius, unsigned &numpanels) {
    int numpanelseach = (numpanels + 9) / 10;
    numpanels = 10 * numpanelseach; // round up

    // Corner points
    Eigen::RowVectorXd x0(2); x0 << outerradius, 0.;
    Eigen::RowVectorXd x1(2); x1 << innerradius*cos(pififths), innerradius*sin(pififths);
    Eigen::RowVectorXd x2(2); x2 << outerradius*cos(twopififths), outerradius*sin(twopififths);
    Eigen::RowVectorXd x3(2); x3 << -innerradius*cos(twopififths), innerradius*sin(twopififths);
    Eigen::RowVectorXd x4(2); x4 << -outerradius*cos(pififths), outerradius*sin(pififths);
    Eigen::RowVectorXd x5(2); x5 << -innerradius, 0.;
    Eigen::RowVectorXd x6(2); x6 << -outerradius*cos(pififths), -outerradius*sin(pififths);
    Eigen::RowVectorXd x7(2); x7 << -innerradius*cos(twopififths), -innerradius*sin(twopififths);
    Eigen::RowVectorXd x8(2); x8 << outerradius*cos(twopififths), -outerradius*sin(twopififths);
    Eigen::RowVectorXd x9(2); x9 << innerradius*cos(pififths), -innerradius*sin(pififths);

    // parametrized line segments forming the edges of the polygon
    ParametrizedLine line0(x0, x1); PanelVector line0panels = line0.split(numpanelseach);
    ParametrizedLine line1(x1, x2); PanelVector line1panels = line1.split(numpanelseach);
    ParametrizedLine line2(x2, x3); PanelVector line2panels = line2.split(numpanelseach);
    ParametrizedLine line3(x3, x4); PanelVector line3panels = line3.split(numpanelseach);
    ParametrizedLine line4(x4, x5); PanelVector line4panels = line4.split(numpanelseach);
    ParametrizedLine line5(x5, x6); PanelVector line5panels = line5.split(numpanelseach);
    ParametrizedLine line6(x6, x7); PanelVector line6panels = line6.split(numpanelseach);
    ParametrizedLine line7(x7, x8); PanelVector line7panels = line7.split(numpanelseach);
    ParametrizedLine line8(x8, x9); PanelVector line8panels = line8.split(numpanelseach);
    ParametrizedLine line9(x9, x0); PanelVector line9panels = line9.split(numpanelseach);

    // storing all the panels in order so that they form a polygon
    PanelVector panels;
    panels.insert(panels.end(), line0panels.begin(), line0panels.end());
    panels.insert(panels.end(), line1panels.begin(), line1panels.end());
    panels.insert(panels.end(), line2panels.begin(), line2panels.end());
    panels.insert(panels.end(), line3panels.begin(), line3panels.end());
    panels.insert(panels.end(), line4panels.begin(), line4panels.end());
    panels.insert(panels.end(), line5panels.begin(), line5panels.end());
    panels.insert(panels.end(), line6panels.begin(), line6panels.end());
    panels.insert(panels.end(), line7panels.begin(), line7panels.end());
    panels.insert(panels.end(), line8panels.begin(), line8panels.end());
    panels.insert(panels.end(), line9panels.begin(), line9panels.end());
    return panels;
}

double normalStar(double innerradius, double outerradius, double x1, double x2) {
    double angle = atan2(x2, x1);
    int wedge_index = floor(angle / pififths);
    double wedge_left = wedge_index * pififths;
    double wedge_right = wedge_left + pififths;

    if (angle - wedge_left < eps) return wedge_left;
    if (wedge_right - angle < eps) return wedge_right;
    if (wedge_index % 2 == 0) // outerradius to left, innerradius to right
        return atan2(outerradius*cos(wedge_left)-innerradius*cos(wedge_right),
                     innerradius*sin(wedge_right)-outerradius*sin(wedge_left));
    // innerradius to left, outerradius to right
    return atan2(innerradius*cos(wedge_left)-outerradius*cos(wedge_right),
                 outerradius*sin(wedge_right)-innerradius*sin(wedge_left));
}

PanelVector buildCShape(double sidelength, double innerradius, unsigned &numpanels) {
    int numpanelseach = (numpanels + 5) / 6;
    numpanels = 6 * numpanelseach; // round up
    
    int numpanelsendcaps = numpanelseach * (innerradius / sidelength);
    if (numpanelsendcaps == 0) numpanelsendcaps = 1;
    
    double halfsidelength = .5 * sidelength;
    // Corner points
    Eigen::RowVectorXd x0(2); x0 << innerradius, halfsidelength - innerradius;
    Eigen::RowVectorXd x1(2); x1 << sidelength + innerradius, halfsidelength - innerradius;
    Eigen::RowVectorXd x2(2); x2 << sidelength + innerradius, halfsidelength + innerradius;
    Eigen::RowVectorXd x3(2); x3 << - innerradius, halfsidelength + innerradius;
    Eigen::RowVectorXd x4(2); x4 << - innerradius, - halfsidelength - innerradius;
    Eigen::RowVectorXd x5(2); x5 << sidelength + innerradius, - halfsidelength - innerradius;
    Eigen::RowVectorXd x6(2); x6 << sidelength + innerradius, - halfsidelength + innerradius;
    Eigen::RowVectorXd x7(2); x7 << innerradius, - halfsidelength + innerradius;

    // parametrized line segments forming the edges of the polygon
    ParametrizedLine line0(x0, x1); PanelVector line0panels = line0.split(numpanelseach);
    ParametrizedLine line1(x1, x2); PanelVector line1panels = line1.split(numpanelsendcaps);
    ParametrizedLine line2(x2, x3); PanelVector line2panels = line2.split(numpanelseach);
    ParametrizedLine line3(x3, x4); PanelVector line3panels = line3.split(numpanelseach);
    ParametrizedLine line4(x4, x5); PanelVector line4panels = line4.split(numpanelseach);
    ParametrizedLine line5(x5, x6); PanelVector line5panels = line5.split(numpanelsendcaps);
    ParametrizedLine line6(x6, x7); PanelVector line6panels = line6.split(numpanelseach);
    ParametrizedLine line7(x7, x0); PanelVector line7panels = line7.split(numpanelseach - 2*numpanelsendcaps);

    // storing all the panels in order so that they form a polygon
    PanelVector panels;
    panels.insert(panels.end(), line0panels.begin(), line0panels.end());
    panels.insert(panels.end(), line1panels.begin(), line1panels.end());
    panels.insert(panels.end(), line2panels.begin(), line2panels.end());
    panels.insert(panels.end(), line3panels.begin(), line3panels.end());
    panels.insert(panels.end(), line4panels.begin(), line4panels.end());
    panels.insert(panels.end(), line5panels.begin(), line5panels.end());
    panels.insert(panels.end(), line6panels.begin(), line6panels.end());
    panels.insert(panels.end(), line7panels.begin(), line7panels.end());
    return panels;
}

double normalCShape(double sidelength, double innerradius, double x1, double x2) {
    double halfsidelength = .5 * sidelength;
    double angleTopL = atan2(x2 - halfsidelength, x1);
    double angleTopR = atan2(x2 - halfsidelength, x1 - sidelength);
    double angleBotL = atan2(x2 + halfsidelength, x1);
    double angleBotR = atan2(x2 + halfsidelength, x1 - sidelength);

    if ( fabs(angleTopR) < piquarters - eps // top endcap
      || fabs(angleBotR) < piquarters - eps // bot endcap
      || ( x1 < innerradius + eps && angleTopL < - piquarters - eps && angleBotL > piquarters + eps )) // right vertical side
        return 0.;
    
    if (( angleTopR > piquarters + eps && angleTopL < threepiquarters - eps ) // top horizontal top side
     || ( x2 < 0. && angleBotR > piquarters + eps && angleBotL < piquarters - eps )) // top horizontal bot side
        return pihalves;
    
    if (( angleBotR < - piquarters - eps && angleBotL > - threepiquarters + eps ) // bot horizontal bot side
     || ( x2 > 0. && angleTopR < - piquarters - eps && angleTopL > - piquarters + eps )) // bot horizontal top side
        return - pihalves;
    
    if ( x1 < 0. && ( angleTopL > threepiquarters + eps || angleTopL < 0. ) && ( angleBotL > 0. || angleBotL < - threepiquarters - eps )) // left vertical side
        return M_PI;
    
    if ( fabs(angleTopR - piquarters) <= eps // top right top endcap
      || fabs(angleBotR - piquarters) <= eps // top right bot endcap
      || ( x1 < innerradius + eps && fabs(angleBotL - piquarters) <= eps )) // bottom right vertical side
        return piquarters;
    
    if ( fabs(angleTopR + piquarters) <= eps // bottom right top endcap
      || fabs(angleBotR + piquarters) <= eps // bottom right bot endcap
      || ( x1 < innerradius + eps && fabs(angleTopL + piquarters) <= eps )) // top right vertical side
        return - piquarters;
    
    if ( fabs(angleTopL - threepiquarters) <= eps ) // top left vertical side
        return threepiquarters;
    
    //if ( fabs(angleBotL + threepiquarters) <= eps ) // bot left vertical side
        return - threepiquarters;
}

PanelVector buildBarbedCShape(double sidelength, double innerradius, unsigned &numpanels) {
    int numpanelseach = (numpanels + 15) / 16;
    numpanels = 16 * numpanelseach; // round up
    
    int numpanelsendcaps = numpanelseach * (innerradius / sidelength);
    if (numpanelsendcaps == 0) numpanelsendcaps = 1;
    
    double halfsidelength = .5 * sidelength;
    // Corner points
    Eigen::RowVectorXd x0(2); x0 << innerradius, halfsidelength - innerradius;
    Eigen::RowVectorXd x1(2); x1 << sidelength - innerradius, halfsidelength - innerradius;
    Eigen::RowVectorXd x2(2); x2 << sidelength - innerradius, innerradius;
    Eigen::RowVectorXd x3(2); x3 << sidelength + innerradius, innerradius;
    Eigen::RowVectorXd x4(2); x4 << sidelength + innerradius, halfsidelength + innerradius;
    Eigen::RowVectorXd x5(2); x5 << - innerradius, halfsidelength + innerradius;
    Eigen::RowVectorXd x6(2); x6 << - innerradius, - halfsidelength - innerradius;
    Eigen::RowVectorXd x7(2); x7 << sidelength + innerradius, - halfsidelength - innerradius;
    Eigen::RowVectorXd x8(2); x8 << sidelength + innerradius, - innerradius;
    Eigen::RowVectorXd x9(2); x9 << sidelength - innerradius, - innerradius;
    Eigen::RowVectorXd xA(2); xA << sidelength - innerradius, - halfsidelength + innerradius;
    Eigen::RowVectorXd xB(2); xB << innerradius, - halfsidelength + innerradius;

    // parametrized line segments forming the edges of the polygon
    ParametrizedLine line0(x0, x1); PanelVector line0panels = line0.split(2 * (numpanelseach - numpanelsendcaps));
    ParametrizedLine line1(x1, x2); PanelVector line1panels = line1.split(numpanelseach - 2 * numpanelsendcaps);
    ParametrizedLine line2(x2, x3); PanelVector line2panels = line2.split(2 * numpanelsendcaps);
    ParametrizedLine line3(x3, x4); PanelVector line3panels = line3.split(numpanelseach);
    ParametrizedLine line4(x4, x5); PanelVector line4panels = line4.split(2 * (numpanelseach + numpanelsendcaps));
    ParametrizedLine line5(x5, x6); PanelVector line5panels = line5.split(2 * (numpanelseach + numpanelsendcaps));
    ParametrizedLine line6(x6, x7); PanelVector line6panels = line6.split(2 * (numpanelseach + numpanelsendcaps));
    ParametrizedLine line7(x7, x8); PanelVector line7panels = line7.split(numpanelseach);
    ParametrizedLine line8(x8, x9); PanelVector line8panels = line8.split(2 * numpanelsendcaps);
    ParametrizedLine line9(x9, xA); PanelVector line9panels = line9.split(numpanelseach - 2 * numpanelsendcaps);
    ParametrizedLine lineA(xA, xB); PanelVector lineApanels = lineA.split(2 * (numpanelseach - numpanelsendcaps));
    ParametrizedLine lineB(xB, x0); PanelVector lineBpanels = lineB.split(2 * (numpanelseach - numpanelsendcaps));

    // storing all the panels in order so that they form a polygon
    PanelVector panels;
    panels.insert(panels.end(), line0panels.begin(), line0panels.end());
    panels.insert(panels.end(), line1panels.begin(), line1panels.end());
    panels.insert(panels.end(), line2panels.begin(), line2panels.end());
    panels.insert(panels.end(), line3panels.begin(), line3panels.end());
    panels.insert(panels.end(), line4panels.begin(), line4panels.end());
    panels.insert(panels.end(), line5panels.begin(), line5panels.end());
    panels.insert(panels.end(), line6panels.begin(), line6panels.end());
    panels.insert(panels.end(), line7panels.begin(), line7panels.end());
    panels.insert(panels.end(), line8panels.begin(), line8panels.end());
    panels.insert(panels.end(), line9panels.begin(), line9panels.end());
    panels.insert(panels.end(), lineApanels.begin(), lineApanels.end());
    panels.insert(panels.end(), lineBpanels.begin(), lineBpanels.end());
    return panels;
}

double normalBarbedCShape(double sidelength, double innerradius, double x1, double x2) {
    double halfsidelength = .5 * sidelength;
    double angleTopL = atan2(x2 - halfsidelength, x1);
    double angleTopR = atan2(x2 - halfsidelength, x1 - sidelength);
    double angleTopER = atan2(x2 - 2. * innerradius, x1 - sidelength);
    double angleBotL = atan2(x2 + halfsidelength, x1);
    double angleBotR = atan2(x2 + halfsidelength, x1 - sidelength);
    double angleBotER = atan2(x2 + 2. * innerradius, x1 - sidelength);

    if ( (angleTopR < piquarters - eps && angleTopER > - piquarters + eps) // right top endcap
      || (angleBotER < piquarters - eps && angleBotR > - piquarters + eps) // right bot endcap
      || ( x1 < innerradius + eps && angleTopL < - piquarters - eps && angleBotL > piquarters + eps )) // right vertical side
        return 0.;
    
    if (( angleTopR > piquarters + eps && angleTopL < threepiquarters - eps ) // top horizontal top side
     || ( x2 < 0. && fabs(angleBotER - pihalves) < piquarters - eps ) // top cap bot side
     || ( x2 < 0. && angleBotR > piquarters + eps && angleBotL < piquarters - eps )) // top horizontal bot side
        return pihalves;
    
    if (( angleBotR < - piquarters - eps && angleBotL > - threepiquarters + eps ) // bot horizontal bot side
     || ( x2 > 0. && fabs(angleTopER + pihalves) < piquarters - eps ) // bot cap top side
     || ( x2 > 0. && angleTopR < - piquarters - eps && angleTopL > - piquarters + eps )) // bot horizontal top side
        return - pihalves;
    
    if (( x1 < 0. && ( angleTopL > threepiquarters + eps || angleTopL < 0. ) && ( angleBotL > 0. || angleBotL < - threepiquarters - eps )) // left vertical side
     || ( x1 > 0. && angleTopR > - threepiquarters + eps && ( angleTopER > threepiquarters + eps || angleTopER < - threepiquarters - eps )) // left cap top side
     || ( x1 > 0. && angleBotR < threepiquarters - eps && ( angleBotER > threepiquarters + eps || angleBotER < - threepiquarters - eps ))) // left cap bot side
        return M_PI;
    
    if ( fabs(angleTopR - piquarters) <= eps // top right top endcap
      || fabs(angleBotER - piquarters) <= eps // top right bot endcap
      || ( x1 < innerradius + eps && fabs(angleBotL - piquarters) <= eps )) // bottom right vertical side
        return piquarters;
    
    if ( fabs(angleTopER + piquarters) <= eps // bottom right top endcap
      || fabs(angleBotR + piquarters) <= eps // bottom right bot endcap
      || ( x1 < innerradius + eps && fabs(angleTopL + piquarters) <= eps )) // top right vertical side
        return - piquarters;
    
    if ( ( x1 > sidelength - innerradius - eps && fabs(angleBotER - threepiquarters) <= eps ) // top left bot endcap
      || ( x1 > sidelength - innerradius - eps && fabs(angleBotER + threepiquarters) <= eps ) // bottom left bot endcap
      || fabs(angleTopL - threepiquarters) <= eps ) // top left vertical side
        return threepiquarters;
    
    //if ( ( x1 > sidelength - innerradius - eps && fabs(angleTopER - threepiquarters) <= eps ) // top left top endcap
    //  || ( x1 > sidelength - innerradius - eps && fabs(angleTopER + threepiquarters) <= eps ) // bottom left top endcap
    //  || fabs(angleTopL - threepiquarters) <= eps ) // bot left vertical side
        return - threepiquarters;
}

PanelVector buildKite(double overhang, double height, unsigned &numpanels) {
    double dt = 2. * M_PI / numpanels;
    PanelVector panels;
    Eigen::RowVectorXd x0(2), x1(2);
    x1 << cos(-dt) + .5 * overhang * (cos(2. * -dt) - 1), height * sin(-dt);
    for (unsigned i = 0; i < numpanels; i++){
        double t = dt * i;
        x0 = x1;
        x1 << cos(t) + .5 * overhang * (cos(2. * t) - 1), height * sin(t);
        ParametrizedLine line(x0, x1);
        PanelVector panel = line.split(1);
        panels.insert(panels.end(), panel.begin(), panel.end());
    }
    return panels;
}

double normalKite(double overhang, double height, double x1, double x2) {
    double theta;
    if (x2 > height)
        theta = pihalves;
    else if (x2 < - height)
        theta = - pihalves;
    else
        theta = asin(x2 / height);

    // angle could be wrong due to left-right symmetry
    if ( height * x1 + overhang * fabs(x2) < 0 )
        theta = M_PI - theta;

    double xdot = - sin(theta) - overhang * sin(2. * theta), ydot = height * cos(theta);
    return atan2(- xdot, ydot);
}

