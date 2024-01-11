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
    numpanels = 4 * ((numpanels + 3) / 4);

    // Corner points
    Eigen::RowVectorXd x1(2); x1 << halfsidelength, -halfsidelength;
    Eigen::RowVectorXd x2(2); x2 << halfsidelength, halfsidelength;
    Eigen::RowVectorXd x3(2); x3 << - halfsidelength, halfsidelength;
    Eigen::RowVectorXd x4(2); x4 << - halfsidelength, - halfsidelength;

    // parametrized line segments forming the edges of the polygon
    ParametrizedLine line1(x1, x2); PanelVector line1panels = line1.split(numpanels/4);
    ParametrizedLine line2(x2, x3); PanelVector line2panels = line2.split(numpanels/4);
    ParametrizedLine line3(x3, x4); PanelVector line3panels = line3.split(numpanels/4);
    ParametrizedLine line4(x4, x1); PanelVector line4panels = line4.split(numpanels/4);

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
    numpanels = 10 * ((numpanels + 9) / 10);

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
    ParametrizedLine line0(x0, x1); PanelVector line0panels = line0.split(numpanels/10);
    ParametrizedLine line1(x1, x2); PanelVector line1panels = line1.split(numpanels/10);
    ParametrizedLine line2(x2, x3); PanelVector line2panels = line2.split(numpanels/10);
    ParametrizedLine line3(x3, x4); PanelVector line3panels = line3.split(numpanels/10);
    ParametrizedLine line4(x4, x5); PanelVector line4panels = line4.split(numpanels/10);
    ParametrizedLine line5(x5, x6); PanelVector line5panels = line5.split(numpanels/10);
    ParametrizedLine line6(x6, x7); PanelVector line6panels = line6.split(numpanels/10);
    ParametrizedLine line7(x7, x8); PanelVector line7panels = line7.split(numpanels/10);
    ParametrizedLine line8(x8, x9); PanelVector line8panels = line8.split(numpanels/10);
    ParametrizedLine line9(x9, x0); PanelVector line9panels = line9.split(numpanels/10);

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
