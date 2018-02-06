#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <QString>

class Parameters
{
public:
    Parameters(int nx, int nt, double range_x, double range_t);

    int get_nx() const;
    int get_nt() const;
    double get_dx() const;
    double get_dt() const;
    double get_alpha() const;

    void set_nx(double nx);
    void set_nt(double nt);
    void set_range_x(double range_x);
    void set_range_t(double range_t);

    QString toQString() const;

private:
    int nx_, nt_;
    double range_x_, range_t_;
    double alpha_;

    void set_alpha();
};

#endif // PARAMETERS_H
