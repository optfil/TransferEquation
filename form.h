#ifndef FORM_H
#define FORM_H

#include <vector>

#include <QComboBox>
#include <QPushButton>
#include <QSlider>
#include <QSpinBox>
#include <QTabWidget>
#include <QTimer>
#include <QWidget>

#include <QtCharts/QtCharts>
QT_CHARTS_USE_NAMESPACE

#include "parameters.h"

constexpr double kRangeX = 10.0;
constexpr double kRangeT = 5.0;
constexpr int kNxMin = 16;
constexpr int kNxMax = 128;
constexpr int kNtMin = 10;
constexpr int kNtMax = 100;

class Form : public QWidget
{
    Q_OBJECT

public:
    Form(QWidget *parent = 0);
    ~Form();

    enum InitialProfile {Gauss, SuperGauss, Rectangle, Step};
    Q_ENUM(InitialProfile)

    enum MethodType {Upwind, Lax, LaxWendroff};
    Q_ENUM(MethodType)

private slots:
    void update_nx_from_slider(int log_n);
    void update_nx(int n);
    void update_nt(int n);
    void updateLabels();
    void initiateState();
    void updateDispersionDiffusionSolution();
    void Solve();
    void Tick();

private:
    QChartView *chartView;
    QLabel *labelInitial;
    QComboBox *comboBoxInitial;
    QLabel *labelSizeX_1, *labelSizeX_2, *labelSizeT_1, *labelSizeT_2, *labelNX_1, *labelNX_2, *labelNT_1, *labelNT_2;
    QLabel *labelSizeX, *labelSizeT;
    QSlider *sliderNX, *sliderNT;
    QSpinBox *spinBoxNX, *spinBoxNT;
    QLabel *labelStepX_1, *labelStepX_2, *labelStepX;
    QLabel *labelStepT_1, *labelStepT_2, *labelStepT;
    QLabel *labelCFL_1, *labelCFL_2, *labelCFL;
    QPushButton *pushButtonSolve;
    QTabWidget *tabWidgetMethods;
    QWidget *upwindWidget, *laxWidget, *laxWendroffWidget;
    QChartView *upwindDispersion, *upwindDissipation, *upwindSolution;
    QChartView *laxDispersion, *laxDissipation, *laxSolution;
    QChartView *laxWendroffDispersion, *laxWendroffDissipation, *laxWendroffSolution;
    QLineSeries *seriesInitial;
    QLineSeries *seriesUpwindIdealDispersion, *seriesLaxIdealDispersion, *seriesLaxWendroffIdealDispersion;
    QLineSeries *seriesUpwindIdealDissipation, *seriesLaxIdealDissipation, *seriesLaxWendroffIdealDissipation;
    QLineSeries *seriesUpwindDispersion, *seriesLaxDispersion, *seriesLaxWendroffDispersion;
    QLineSeries *seriesUpwindDissipation, *seriesLaxDissipation, *seriesLaxWendroffDissipation;

    QTimer *timer;

    Parameters *param;
    MethodType method_;
    std::vector<double> state_;
    std::vector<double> tmp_state_;
    double t_cur_;

    void showState();
};

#endif // FORM_H
