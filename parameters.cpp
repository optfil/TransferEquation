#include "parameters.h"

Parameters::Parameters(int nx, int nt, double range_x, double range_t)
    : nx_(nx), nt_(nt), range_x_(range_x), range_t_(range_t)
{
    set_alpha();
}

int Parameters::get_nx() const
{
    return nx_;
}

int Parameters::get_nt() const
{
    return nt_;
}

double Parameters::get_dx() const
{
    return range_x_ / (nx_-1);
}

double Parameters::get_dt() const
{
    return range_t_ / nt_;
}

double Parameters::get_alpha() const
{
    return alpha_;
}

void Parameters::set_nx(double nx)
{
    nx_ = nx;
    set_alpha();
}

void Parameters::set_nt(double nt)
{
    nt_ = nt;
    set_alpha();
}

void Parameters::set_range_x(double range_x)
{
    range_x_ = range_x;
    set_alpha();
}

void Parameters::set_range_t(double range_t)
{
    range_t_ = range_t;
    set_alpha();
}

void Parameters::set_alpha()
{
    alpha_ = get_dt() / get_dx();
}

QString Parameters::toQString() const
{
    return QString::number(nx_) + " " + QString::number(nt_);
}
