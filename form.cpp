#include "form.h"

#include <QLayout>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <utility>

#include "fftw3.h"


static double initial(double x, Form::InitialProfile profile)
{
    switch (profile)
    {
    case Form::Gauss:
        return std::exp(-std::pow(x - kRangeX/4.0, 2.0));
    case Form::SuperGauss:
        return std::exp(-std::pow(x - kRangeX/4.0, 8.0));
    case Form::Rectangle:
        return (x > kRangeX/8.0 && x < kRangeX * 3.0 / 8.0) ? 1.0 : 0.0;
    case Form::Step:
        return (x < kRangeX/4.0) ? 0.0 : 1.0;
    default:
        return 0;
    }
}

static std::pair<double, double> dispersion_diffusion(double q_N, double alpha, Form::MethodType type)
{
    std::complex<double> lambda;
    double kappa = 2.0*M_PI*q_N;
    switch (type)
    {
    case Form::Upwind:
        lambda = 1.0 - alpha * (1.0 - std::exp(std::complex<double>(0.0, kappa)));
        break;
    case Form::Lax:
        lambda = std::complex<double>(std::cos(kappa), alpha * std::sin(kappa));
        break;
    case Form::LaxWendroff:
        lambda = std::complex<double>(1.0 - alpha*alpha * (1.0 - std::cos(kappa)), alpha * std::sin(kappa));
        break;
    default:
        lambda = 1.0;
        break;
    }

    lambda = std::log(lambda);

    return std::make_pair(std::imag(lambda), -std::real(lambda));
}

static void setGrid(QValueAxis* ax)
{
    ax->setGridLineVisible(true);
    QPen pen = ax->gridLinePen();
    pen.setWidth(2);
    pen.setColor(Qt::gray);
    ax->setGridLinePen(pen);
}

Form::Form(QWidget *parent)
    : QWidget(parent), param(nullptr), t_cur_(0.0)
{
    timer = new QTimer();
    timer->setInterval(30);

    seriesInitial = new QLineSeries();
    seriesInitial->setColor(Qt::blue);
    seriesInitial->setPen(QPen(seriesInitial->pen().brush(), 3));
    seriesUpwindIdealDispersion = new QLineSeries();
    seriesUpwindIdealDispersion->setColor(Qt::blue);
    seriesUpwindIdealDispersion->setPen(QPen(seriesUpwindIdealDispersion->pen().brush(), 3));
    seriesLaxIdealDispersion = new QLineSeries();
    seriesLaxIdealDispersion->setColor(Qt::blue);
    seriesLaxIdealDispersion->setPen(QPen(seriesLaxIdealDispersion->pen().brush(), 3));
    seriesLaxWendroffIdealDispersion = new QLineSeries();
    seriesLaxWendroffIdealDispersion->setColor(Qt::blue);
    seriesLaxWendroffIdealDispersion->setPen(QPen(seriesLaxWendroffIdealDispersion->pen().brush(), 3));
    seriesUpwindIdealDissipation = new QLineSeries();
    seriesUpwindIdealDissipation->setColor(Qt::blue);
    seriesUpwindIdealDissipation->setPen(QPen(seriesUpwindIdealDissipation->pen().brush(), 3));
    seriesLaxIdealDissipation = new QLineSeries();
    seriesLaxIdealDissipation->setColor(Qt::blue);
    seriesLaxIdealDissipation->setPen(QPen(seriesLaxIdealDissipation->pen().brush(), 3));
    seriesLaxWendroffIdealDissipation = new QLineSeries();
    seriesLaxWendroffIdealDissipation->setColor(Qt::blue);
    seriesLaxWendroffIdealDissipation->setPen(QPen(seriesLaxWendroffIdealDissipation->pen().brush(), 3));
    seriesUpwindDispersion = new QLineSeries();
    seriesUpwindDispersion->setColor(Qt::red);
    seriesUpwindDispersion->setPen(QPen(seriesUpwindDispersion->pen().brush(), 3));
    seriesLaxDispersion = new QLineSeries();
    seriesLaxDispersion->setColor(Qt::red);
    seriesLaxDispersion->setPen(QPen(seriesLaxDispersion->pen().brush(), 3));
    seriesLaxWendroffDispersion = new QLineSeries();
    seriesLaxWendroffDispersion->setColor(Qt::red);
    seriesLaxWendroffDispersion->setPen(QPen(seriesLaxWendroffDispersion->pen().brush(), 3));
    seriesUpwindDissipation = new QLineSeries();
    seriesUpwindDissipation->setColor(Qt::red);
    seriesUpwindDissipation->setPen(QPen(seriesUpwindDissipation->pen().brush(), 3));
    seriesLaxDissipation = new QLineSeries();
    seriesLaxDissipation->setColor(Qt::red);
    seriesLaxDissipation->setPen(QPen(seriesLaxDissipation->pen().brush(), 3));
    seriesLaxWendroffDissipation = new QLineSeries();
    seriesLaxWendroffDissipation->setColor(Qt::red);
    seriesLaxWendroffDissipation->setPen(QPen(seriesLaxWendroffDissipation->pen().brush(), 3));
    spectreUpwindDispersion = new QBarSeries();
    spectreUpwindDissipation = new QBarSeries();
    spectreLaxDispersion = new QBarSeries();
    spectreLaxDissipation = new QBarSeries();
    spectreLaxWendroffDispersion = new QBarSeries();
    spectreLaxWendroffDissipation = new QBarSeries();

    seriesUpwindIdealDissipation->append(QList<QPointF>() << QPointF(0.0, 0.0) << QPointF(0.5, 0.0));
    seriesLaxIdealDissipation->append(QList<QPointF>() << QPointF(0.0, 0.0) << QPointF(0.5, 0.0));
    seriesLaxWendroffIdealDissipation->append(QList<QPointF>() << QPointF(0.0, 0.0) << QPointF(0.5, 0.0));

    QChart *chartInitial = new QChart();
    chartInitial->addSeries(seriesInitial);
    chartInitial->setTitle(tr("Initial profile"));
    chartInitial->legend()->hide();

    QValueAxis *axisXInitial = new QValueAxis;
    axisXInitial->setLineVisible(false);
    setGrid(axisXInitial);
    axisXInitial->setLabelsVisible(false);
    axisXInitial->setRange(0.0, kRangeX/2);
    chartInitial->addAxis(axisXInitial, Qt::AlignBottom);
    seriesInitial->attachAxis(axisXInitial);
    QValueAxis *axisYInitial = new QValueAxis;
    axisYInitial->setLineVisible(false);
    setGrid(axisYInitial);
    axisYInitial->setLabelsVisible(false);
    axisYInitial->setRange(0.0, 1.003);
    chartInitial->addAxis(axisYInitial, Qt::AlignLeft);
    seriesInitial->attachAxis(axisYInitial);

    chartView = new QChartView(chartInitial);
    chartView->setRenderHint(QPainter::Antialiasing);

    labelInitial = new QLabel(tr("Pulse form"));
    comboBoxInitial = new QComboBox();
    comboBoxInitial->addItem(tr("Gauss"), QVariant(Gauss));
    comboBoxInitial->addItem(tr("SuperGauss"), QVariant(SuperGauss));
    comboBoxInitial->addItem(tr("Rectangle"), QVariant(Rectangle));
    comboBoxInitial->addItem(tr("Step"), QVariant(Step));

    labelSizeX_1 = new QLabel(tr("Grid size"));
    labelSizeX_2 = new QLabel(tr(" L = "));
    labelSizeX_2->setAlignment(Qt::AlignRight);
    labelSizeT_1 = new QLabel(tr("Integration time"));
    labelSizeT_2 = new QLabel(tr(" T = "));
    labelSizeT_2->setAlignment(Qt::AlignRight);
    labelNX_1 = new QLabel(tr("Number of spatial points"));
    labelNX_2 = new QLabel(tr("NX = "));
    labelNX_2->setAlignment(Qt::AlignRight);
    labelNT_1 = new QLabel(tr("Number of temporal points"));
    labelNT_2 = new QLabel(tr("NT = "));
    labelNT_2->setAlignment(Qt::AlignRight);

    labelSizeX = new QLabel(QString::number(kRangeX, 'f', 1));
    labelSizeT = new QLabel(QString::number(kRangeT, 'f', 1));

    sliderNX = new QSlider(Qt::Horizontal);
    sliderNX->setRange(1, 4);
    sliderNX->setSingleStep(1);
    sliderNX->setTickInterval(1);
    sliderNX->setTickPosition(QSlider::TicksBelow);
    sliderNX->setValue(1);

    sliderNT = new QSlider(Qt::Horizontal);
    sliderNT->setRange(kNtMin, kNtMax);
    sliderNT->setTickInterval(10);
    sliderNT->setTickPosition(QSlider::TicksBelow);
    sliderNT->setValue(kNtMin);

    spinBoxNX = new QSpinBox();
    spinBoxNX->setMinimum(0);
    spinBoxNX->setMaximum(kNxMax);
    spinBoxNX->setValue(kNxMin);
    spinBoxNX->setSingleStep(kNxMin);

    spinBoxNT = new QSpinBox();
    spinBoxNT->setMinimum(kNtMin);
    spinBoxNT->setMaximum(kNtMax);
    spinBoxNT->setValue(kNtMin);
    spinBoxNT->setSingleStep(1);

    labelStepX_1 = new QLabel(tr("Spatial step"));
    labelStepX_2 = new QLabel(tr("dx = "));
    labelStepX_2->setAlignment(Qt::AlignRight);
    labelStepX = new QLabel();

    labelStepT_1 = new QLabel(tr("Time step"));
    labelStepT_2 = new QLabel(tr("dt = "));
    labelStepT_2->setAlignment(Qt::AlignRight);
    labelStepT = new QLabel();

    labelCFL_1 = new QLabel(tr("CFL number"));
    labelCFL_2 = new QLabel(tr("α = "));
    labelCFL_2->setAlignment(Qt::AlignRight);
    labelCFL = new QLabel();

    pushButtonSolve = new QPushButton(tr("Start"));

    upwindWidget = new QWidget();

    QChart *upwindDispersionChart = new QChart();
    upwindDispersionChart->addSeries(seriesUpwindIdealDispersion);
    upwindDispersionChart->addSeries(seriesUpwindDispersion);
    upwindDispersionChart->addSeries(spectreUpwindDispersion);
    upwindDispersionChart->setTitle(tr("Dispersion error"));
    upwindDispersionChart->legend()->hide();

    QValueAxis *axisXUpwindDispersion = new QValueAxis;
    axisXUpwindDispersion->setLineVisible(false);
    setGrid(axisXUpwindDispersion);
    axisXUpwindDispersion->setLabelsVisible(false);
    axisXUpwindDispersion->setRange(0.0, 0.5);
    upwindDispersionChart->addAxis(axisXUpwindDispersion, Qt::AlignBottom);
    seriesUpwindIdealDispersion->attachAxis(axisXUpwindDispersion);
    seriesUpwindDispersion->attachAxis(axisXUpwindDispersion);
    spectreUpwindDispersion->attachAxis(axisXUpwindDispersion);
    QValueAxis *axisYUpwindDispersion = new QValueAxis;
    axisYUpwindDispersion->setLineVisible(false);
    setGrid(axisYUpwindDispersion);
    axisYUpwindDispersion->setLabelsVisible(false);
    axisYUpwindDispersion->setRange(0.0, 2.0);
    upwindDispersionChart->addAxis(axisYUpwindDispersion, Qt::AlignLeft);
    seriesUpwindIdealDispersion->attachAxis(axisYUpwindDispersion);
    seriesUpwindDispersion->attachAxis(axisYUpwindDispersion);
    spectreUpwindDispersion->attachAxis(axisYUpwindDispersion);

    upwindDispersion = new QChartView();
    upwindDispersion->setRenderHint(QPainter::Antialiasing);
    upwindDispersion->setChart(upwindDispersionChart);

    QChart *upwindDissipationChart = new QChart();
    upwindDissipationChart->addSeries(seriesUpwindIdealDissipation);
    upwindDissipationChart->addSeries(seriesUpwindDissipation);
    upwindDissipationChart->addSeries(spectreUpwindDissipation);
    upwindDissipationChart->setTitle(tr("Dissipation error"));
    upwindDissipationChart->legend()->hide();

    QValueAxis *axisXUpwindDissipation = new QValueAxis;
    axisXUpwindDissipation->setLineVisible(false);
    setGrid(axisXUpwindDissipation);
    axisXUpwindDissipation->setLabelsVisible(false);
    axisXUpwindDissipation->setRange(0.0, 0.5);
    upwindDissipationChart->addAxis(axisXUpwindDissipation, Qt::AlignBottom);
    seriesUpwindIdealDissipation->attachAxis(axisXUpwindDissipation);
    seriesUpwindDissipation->attachAxis(axisXUpwindDissipation);
    spectreUpwindDissipation->attachAxis(axisXUpwindDissipation);
    QValueAxis *axisYUpwindDissipation = new QValueAxis;
    axisYUpwindDissipation->setLineVisible(false);
    setGrid(axisYUpwindDissipation);
    axisYUpwindDissipation->setLabelsVisible(false);
    axisYUpwindDissipation->setRange(-3.0, 3.0);
    upwindDissipationChart->addAxis(axisYUpwindDissipation, Qt::AlignLeft);
    seriesUpwindIdealDissipation->attachAxis(axisYUpwindDissipation);
    seriesUpwindDissipation->attachAxis(axisYUpwindDissipation);
    spectreUpwindDissipation->attachAxis(axisYUpwindDissipation);

    upwindDissipation = new QChartView();
    upwindDissipation->setRenderHint(QPainter::Antialiasing);
    upwindDissipation->setChart(upwindDissipationChart);

    QChart *upwindSolutionChart = new QChart();
    upwindSolutionChart->setTitle(tr("Solution"));
    upwindSolutionChart->legend()->hide();
    QValueAxis *axisXUpwindSolution = new QValueAxis;
    axisXUpwindSolution->setLineVisible(false);
    setGrid(axisXUpwindSolution);
    axisXUpwindSolution->setLabelsVisible(false);
    axisXUpwindSolution->setRange(0.0, kRangeX);
    upwindSolutionChart->addAxis(axisXUpwindSolution, Qt::AlignBottom);
    QValueAxis *axisYUpwindSolution = new QValueAxis;
    axisYUpwindSolution->setLineVisible(false);
    setGrid(axisYUpwindSolution);
    axisYUpwindSolution->setLabelsVisible(false);
    axisYUpwindSolution->setRange(-0.5, 1.5);
    upwindSolutionChart->addAxis(axisYUpwindSolution, Qt::AlignLeft);

    upwindSolution = new QChartView();
    upwindSolution->setRenderHint(QPainter::Antialiasing);
    upwindSolution->setChart(upwindSolutionChart);

    QVBoxLayout *upwindLeft = new QVBoxLayout();
    upwindLeft->addWidget(upwindDispersion);
    upwindLeft->addWidget(upwindDissipation);
    QHBoxLayout *upwindMain = new QHBoxLayout();
    upwindMain->addLayout(upwindLeft);
    upwindMain->addWidget(upwindSolution);
    upwindWidget->setLayout(upwindMain);

    laxWidget = new QWidget();

    QChart *laxDispersionChart = new QChart();
    laxDispersionChart->addSeries(seriesLaxIdealDispersion);
    laxDispersionChart->addSeries(seriesLaxDispersion);
    laxDispersionChart->addSeries(spectreLaxDispersion);
    laxDispersionChart->setTitle(tr("Dispersion error"));
    laxDispersionChart->legend()->hide();

    QValueAxis *axisXLaxDispersion = new QValueAxis;
    axisXLaxDispersion->setLineVisible(false);
    setGrid(axisXLaxDispersion);
    axisXLaxDispersion->setLabelsVisible(false);
    axisXLaxDispersion->setRange(0.0, 0.5);
    laxDispersionChart->addAxis(axisXLaxDispersion, Qt::AlignBottom);
    seriesLaxIdealDispersion->attachAxis(axisXLaxDispersion);
    seriesLaxDispersion->attachAxis(axisXLaxDispersion);
    spectreLaxDispersion->attachAxis(axisXLaxDispersion);
    QValueAxis *axisYLaxDispersion = new QValueAxis;
    axisYLaxDispersion->setLineVisible(false);
    setGrid(axisYLaxDispersion);
    axisYLaxDispersion->setLabelsVisible(false);
    axisYLaxDispersion->setRange(0.0, 2.0);
    laxDispersionChart->addAxis(axisYLaxDispersion, Qt::AlignLeft);
    seriesLaxIdealDispersion->attachAxis(axisYLaxDispersion);
    seriesLaxDispersion->attachAxis(axisYLaxDispersion);
    spectreLaxDispersion->attachAxis(axisYLaxDispersion);

    laxDispersion = new QChartView();
    laxDispersion->setRenderHint(QPainter::Antialiasing);
    laxDispersion->setChart(laxDispersionChart);

    QChart *laxDissipationChart = new QChart();
    laxDissipationChart->addSeries(seriesLaxIdealDissipation);
    laxDissipationChart->addSeries(seriesLaxDissipation);
    laxDissipationChart->addSeries(spectreLaxDissipation);
    laxDissipationChart->setTitle(tr("Dissipation error"));
    laxDissipationChart->legend()->hide();

    QValueAxis *axisXLaxDissipation = new QValueAxis;
    axisXLaxDissipation->setLineVisible(false);
    setGrid(axisXLaxDissipation);
    axisXLaxDissipation->setLabelsVisible(false);
    axisXLaxDissipation->setRange(0.0, 0.5);
    laxDissipationChart->addAxis(axisXLaxDissipation, Qt::AlignBottom);
    seriesLaxIdealDissipation->attachAxis(axisXLaxDissipation);
    seriesLaxDissipation->attachAxis(axisXLaxDissipation);
    spectreLaxDissipation->attachAxis(axisXLaxDissipation);
    QValueAxis *axisYLaxDissipation = new QValueAxis;
    axisYLaxDissipation->setLineVisible(false);
    setGrid(axisYLaxDissipation);
    axisYLaxDissipation->setLabelsVisible(false);
    axisYLaxDissipation->setRange(-3.0, 3.0);
    laxDissipationChart->addAxis(axisYLaxDissipation, Qt::AlignLeft);
    seriesLaxIdealDissipation->attachAxis(axisYLaxDissipation);
    seriesLaxDissipation->attachAxis(axisYLaxDissipation);
    spectreLaxDissipation->attachAxis(axisYLaxDissipation);

    laxDissipation = new QChartView();
    laxDissipation->setRenderHint(QPainter::Antialiasing);
    laxDissipation->setChart(laxDissipationChart);

    QChart *laxSolutionChart = new QChart();
    laxSolutionChart->setTitle(tr("Solution"));
    laxSolutionChart->legend()->hide();
    QValueAxis *axisXLaxSolution = new QValueAxis;
    axisXLaxSolution->setLineVisible(false);
    setGrid(axisXLaxSolution);
    axisXLaxSolution->setLabelsVisible(false);
    axisXLaxSolution->setRange(0.0, kRangeX);
    laxSolutionChart->addAxis(axisXLaxSolution, Qt::AlignBottom);
    QValueAxis *axisYLaxSolution = new QValueAxis;
    axisYLaxSolution->setLineVisible(false);
    setGrid(axisYLaxSolution);
    axisYLaxSolution->setLabelsVisible(false);
    axisYLaxSolution->setRange(-0.5, 1.5);
    laxSolutionChart->addAxis(axisYLaxSolution, Qt::AlignLeft);

    laxSolution = new QChartView();
    laxSolution->setRenderHint(QPainter::Antialiasing);
    laxSolution->setChart(laxSolutionChart);

    QVBoxLayout *laxLeft = new QVBoxLayout();
    laxLeft->addWidget(laxDispersion);
    laxLeft->addWidget(laxDissipation);
    QHBoxLayout *laxMain = new QHBoxLayout();
    laxMain->addLayout(laxLeft);
    laxMain->addWidget(laxSolution);
    laxWidget->setLayout(laxMain);

    laxWendroffWidget = new QWidget();

    QChart *laxWendroffDispersionChart = new QChart();
    laxWendroffDispersionChart->addSeries(seriesLaxWendroffIdealDispersion);
    laxWendroffDispersionChart->addSeries(seriesLaxWendroffDispersion);
    laxWendroffDispersionChart->addSeries(spectreLaxWendroffDispersion);
    laxWendroffDispersionChart->setTitle(tr("Dispersion error"));
    laxWendroffDispersionChart->legend()->hide();

    QValueAxis *axisXLaxWendroffDispersion = new QValueAxis;
    axisXLaxWendroffDispersion->setLineVisible(false);
    setGrid(axisXLaxWendroffDispersion);
    axisXLaxWendroffDispersion->setLabelsVisible(false);
    axisXLaxWendroffDispersion->setRange(0.0, 0.5);
    laxWendroffDispersionChart->addAxis(axisXLaxWendroffDispersion, Qt::AlignBottom);
    seriesLaxWendroffIdealDispersion->attachAxis(axisXLaxWendroffDispersion);
    seriesLaxWendroffDispersion->attachAxis(axisXLaxWendroffDispersion);
    spectreLaxWendroffDispersion->attachAxis(axisXLaxWendroffDispersion);
    QValueAxis *axisYLaxWendroffDispersion = new QValueAxis;
    axisYLaxWendroffDispersion->setLineVisible(false);
    setGrid(axisYLaxWendroffDispersion);
    axisYLaxWendroffDispersion->setLabelsVisible(false);
    axisYLaxWendroffDispersion->setRange(0.0, 2.0);
    laxWendroffDispersionChart->addAxis(axisYLaxWendroffDispersion, Qt::AlignLeft);
    seriesLaxWendroffIdealDispersion->attachAxis(axisYLaxWendroffDispersion);
    seriesLaxWendroffDispersion->attachAxis(axisYLaxWendroffDispersion);
    spectreLaxWendroffDispersion->attachAxis(axisYLaxWendroffDispersion);

    laxWendroffDispersion = new QChartView();
    laxWendroffDispersion->setRenderHint(QPainter::Antialiasing);
    laxWendroffDispersion->setChart(laxWendroffDispersionChart);

    QChart *laxWendroffDissipationChart = new QChart();
    laxWendroffDissipationChart->addSeries(seriesLaxWendroffIdealDissipation);
    laxWendroffDissipationChart->addSeries(seriesLaxWendroffDissipation);
    laxWendroffDissipationChart->addSeries(spectreLaxWendroffDissipation);
    laxWendroffDissipationChart->setTitle(tr("Dissipation error"));
    laxWendroffDissipationChart->legend()->hide();

    QValueAxis *axisXLaxWendroffDissipation = new QValueAxis;
    axisXLaxWendroffDissipation->setLineVisible(false);
    setGrid(axisXLaxWendroffDissipation);
    axisXLaxWendroffDissipation->setLabelsVisible(false);
    axisXLaxWendroffDissipation->setRange(0.0, 0.5);
    laxWendroffDissipationChart->addAxis(axisXLaxWendroffDissipation, Qt::AlignBottom);
    seriesLaxWendroffIdealDissipation->attachAxis(axisXLaxWendroffDissipation);
    seriesLaxWendroffDissipation->attachAxis(axisXLaxWendroffDissipation);
    spectreLaxWendroffDissipation->attachAxis(axisXLaxWendroffDissipation);
    QValueAxis *axisYLaxWendroffDissipation = new QValueAxis;
    axisYLaxWendroffDissipation->setLineVisible(false);
    setGrid(axisYLaxWendroffDissipation);
    axisYLaxWendroffDissipation->setLabelsVisible(false);
    axisYLaxWendroffDissipation->setRange(-3.0, 3.0);
    laxWendroffDissipationChart->addAxis(axisYLaxWendroffDissipation, Qt::AlignLeft);
    seriesLaxWendroffIdealDissipation->attachAxis(axisYLaxWendroffDissipation);
    seriesLaxWendroffDissipation->attachAxis(axisYLaxWendroffDissipation);
    spectreUpwindDissipation->attachAxis(axisYLaxWendroffDissipation);

    laxWendroffDissipation = new QChartView();
    laxWendroffDissipation->setRenderHint(QPainter::Antialiasing);
    laxWendroffDissipation->setChart(laxWendroffDissipationChart);

    QChart *laxWendroffSolutionChart = new QChart();
    laxWendroffSolutionChart->setTitle(tr("Solution"));
    laxWendroffSolutionChart->legend()->hide();
    QValueAxis *axisXLaxWendroffSolution = new QValueAxis;
    axisXLaxWendroffSolution->setLineVisible(false);
    setGrid(axisXLaxWendroffSolution);
    axisXLaxWendroffSolution->setLabelsVisible(false);
    axisXLaxWendroffSolution->setRange(0.0, kRangeX);
    laxWendroffSolutionChart->addAxis(axisXLaxWendroffSolution, Qt::AlignBottom);
    QValueAxis *axisYLaxWendroffSolution = new QValueAxis;
    axisYLaxWendroffSolution->setLineVisible(false);
    setGrid(axisYLaxWendroffSolution);
    axisYLaxWendroffSolution->setLabelsVisible(false);
    axisYLaxWendroffSolution->setRange(-0.5, 1.5);
    laxWendroffSolutionChart->addAxis(axisYLaxWendroffSolution, Qt::AlignLeft);

    laxWendroffSolution = new QChartView();
    laxWendroffSolution->setRenderHint(QPainter::Antialiasing);
    laxWendroffSolution->setChart(laxWendroffSolutionChart);

    QVBoxLayout *laxWendroffLeft = new QVBoxLayout();
    laxWendroffLeft->addWidget(laxWendroffDispersion);
    laxWendroffLeft->addWidget(laxWendroffDissipation);
    QHBoxLayout *laxWendroffMain = new QHBoxLayout();
    laxWendroffMain->addLayout(laxWendroffLeft);
    laxWendroffMain->addWidget(laxWendroffSolution);
    laxWendroffWidget->setLayout(laxWendroffMain);

    tabWidgetMethods = new QTabWidget();
    tabWidgetMethods->addTab(upwindWidget, tr("Upwind"));
    tabWidgetMethods->addTab(laxWidget, tr("Lax-Friedrichs"));
    tabWidgetMethods->addTab(laxWendroffWidget, tr("Lax-Wendroff"));

    QGridLayout *layoutNxNt = new QGridLayout();
    layoutNxNt->addWidget(labelInitial, 0, 0, 1, 1);
    layoutNxNt->addWidget(comboBoxInitial, 0, 1, 1, 3);
    layoutNxNt->addWidget(labelSizeX_1, 1, 0, 1, 1);
    layoutNxNt->addWidget(labelSizeX_2, 1, 1, 1, 1);
    layoutNxNt->addWidget(labelSizeX, 1, 2, 1, 1);
    layoutNxNt->addWidget(labelSizeT_1, 2, 0, 1, 1);
    layoutNxNt->addWidget(labelSizeT_2, 2, 1, 1, 1);
    layoutNxNt->addWidget(labelSizeT, 2, 2, 1, 1);
    layoutNxNt->addWidget(labelNX_1, 3, 0, 1, 1, Qt::AlignBaseline);
    layoutNxNt->addWidget(labelNX_2, 3, 1, 1, 1, Qt::AlignBaseline);
    layoutNxNt->addWidget(spinBoxNX, 3, 2, 1, 1);
    layoutNxNt->addWidget(sliderNX, 3, 3, 1, 1);
    layoutNxNt->addWidget(labelNT_1, 4, 0, 1, 1, Qt::AlignBaseline);
    layoutNxNt->addWidget(labelNT_2, 4, 1, 1, 1, Qt::AlignBaseline);
    layoutNxNt->addWidget(spinBoxNT, 4, 2, 1, 1);
    layoutNxNt->addWidget(sliderNT, 4, 3, 1, 1);
    layoutNxNt->addWidget(labelStepX_1, 5, 0, 1, 1);
    layoutNxNt->addWidget(labelStepX_2, 5, 1, 1, 1);
    layoutNxNt->addWidget(labelStepX, 5, 2, 1, 1);
    layoutNxNt->addWidget(labelStepT_1, 6, 0, 1, 1);
    layoutNxNt->addWidget(labelStepT_2, 6, 1, 1, 1);
    layoutNxNt->addWidget(labelStepT, 6, 2, 1, 1);
    layoutNxNt->addWidget(labelCFL_1, 7, 0, 1, 1);
    layoutNxNt->addWidget(labelCFL_2, 7, 1, 1, 1);
    layoutNxNt->addWidget(labelCFL, 7, 2, 1, 1);
    layoutNxNt->addWidget(pushButtonSolve, 6, 3, 2, 1);

    QVBoxLayout *layoutParam = new QVBoxLayout;
    layoutParam->addWidget(chartView);
    layoutParam->addLayout(layoutNxNt);

    QHBoxLayout *layoutMain = new QHBoxLayout();
    layoutMain->addLayout(layoutParam);
    layoutMain->addWidget(tabWidgetMethods);

    setLayout(layoutMain);

    connect(comboBoxInitial, SIGNAL(currentIndexChanged(int)), this, SLOT(initiateState()));
    connect(sliderNX, SIGNAL(valueChanged(int)), this, SLOT(update_nx_from_slider(int)));
    connect(sliderNT, SIGNAL(valueChanged(int)), this, SLOT(update_nt(int)));
    connect(spinBoxNX, SIGNAL(valueChanged(int)), this, SLOT(update_nx(int)));
    connect(spinBoxNT, SIGNAL(valueChanged(int)), this, SLOT(update_nt(int)));
    connect(tabWidgetMethods, SIGNAL(currentChanged(int)), this, SLOT(updateDispersionDiffusionSolution()));
    connect(pushButtonSolve, SIGNAL(clicked(bool)), this, SLOT(Solve()));
    connect(timer, SIGNAL(timeout()), this, SLOT(Tick()));

    initiateState();
}

Form::~Form()
{
    delete param;
}

void Form::update_nx_from_slider(int n)
{
    int new_nx = static_cast<int>(std::round(std::pow(2.0, n+3)));

    spinBoxNX->blockSignals(true);
    sliderNX->blockSignals(true);

    spinBoxNX->setValue(new_nx);
    spinBoxNX->setSingleStep(new_nx);

    spinBoxNX->blockSignals(false);
    sliderNX->blockSignals(false);

    initiateState();
}

void Form::update_nx(int n)
{
    int old_nx = param->get_nx();
    if (n == 0)
    {
        n = std::max(old_nx/2, kNxMin);
    }

    int new_nx_log = static_cast<int>(std::round(std::log2(static_cast<double>(n))));
    int new_nx = static_cast<int>(std::round(std::pow(2.0, new_nx_log)));

    spinBoxNX->blockSignals(true);
    sliderNX->blockSignals(true);

    sliderNX->setValue(new_nx_log-3);
    spinBoxNX->setValue(new_nx);
    spinBoxNX->setSingleStep(new_nx);

    spinBoxNX->blockSignals(false);
    sliderNX->blockSignals(false);

    initiateState();
}

void Form::update_nt(int n)
{
    spinBoxNT->blockSignals(true);
    sliderNT->blockSignals(true);

    sliderNT->setValue(n);
    spinBoxNT->setValue(n);

    spinBoxNT->blockSignals(false);
    sliderNT->blockSignals(false);

    initiateState();
}

void Form::updateLabels()
{
    labelStepX->setText(QString::number(param->get_dx(), 'f', 3));
    labelStepT->setText(QString::number(param->get_dt(), 'f', 3));
    labelCFL->setText(QString::number(param->get_alpha(), 'f', 3));
}

void Form::updateDispersionDiffusionSolution()
{
    method_ = static_cast<MethodType>(tabWidgetMethods->currentIndex());

    QList<QPointF> upwind_disp_data, upwind_diff_data, lax_disp_data, lax_diff_data, lax_wendroff_disp_data, lax_wendroff_diff_data;
    std::pair<double, double> coeffs;
    double ideal_disp_max = 2.0*M_PI*param->get_alpha() * 0.5;
    for (int i = 0; i < param->get_nx()/2+1; ++i)
    {
        double xi = static_cast<double>(i) / param->get_nx();

        coeffs = dispersion_diffusion(xi, param->get_alpha(), Form::Upwind);
        upwind_disp_data.append(QPointF(xi, coeffs.first));
        upwind_diff_data.append(QPointF(xi, coeffs.second));

        coeffs = dispersion_diffusion(xi, param->get_alpha(), Form::Lax);
        lax_disp_data.append(QPointF(xi, coeffs.first));
        lax_diff_data.append(QPointF(xi, coeffs.second));

        coeffs = dispersion_diffusion(xi, param->get_alpha(), Form::LaxWendroff);
        lax_wendroff_disp_data.append(QPointF(xi, coeffs.first));
        lax_wendroff_diff_data.append(QPointF(xi, coeffs.second));
    }

    for (int i = 0; i < param->get_nx()/2+1; ++i)
    {
        upwind_disp_data[i].setY(upwind_disp_data[i].y() / ideal_disp_max);
        lax_disp_data[i].setY(lax_disp_data[i].y() / ideal_disp_max);
        lax_wendroff_disp_data[i].setY(lax_wendroff_disp_data[i].y() / ideal_disp_max);
    }

    seriesUpwindIdealDispersion->clear();
    seriesUpwindIdealDispersion->append(QList<QPointF>() << QPointF(0.0, 0.0) << QPointF(0.5, ideal_disp_max / ideal_disp_max));
    seriesUpwindDispersion->clear();
    seriesUpwindDispersion->append(upwind_disp_data);
    seriesUpwindDissipation->clear();
    seriesUpwindDissipation->append(upwind_diff_data);

    seriesLaxIdealDispersion->clear();
    seriesLaxIdealDispersion->append(QList<QPointF>() << QPointF(0.0, 0.0) << QPointF(0.5, ideal_disp_max / ideal_disp_max));
    seriesLaxDispersion->clear();
    seriesLaxDispersion->append(lax_disp_data);
    seriesLaxDissipation->clear();
    seriesLaxDissipation->append(lax_diff_data);

    seriesLaxWendroffIdealDispersion->clear();
    seriesLaxWendroffIdealDispersion->append(QList<QPointF>() << QPointF(0.0, 0.0) << QPointF(0.5, ideal_disp_max / ideal_disp_max));
    seriesLaxWendroffDispersion->clear();
    seriesLaxWendroffDispersion->append(lax_wendroff_disp_data);
    seriesLaxWendroffDissipation->clear();
    seriesLaxWendroffDissipation->append(lax_wendroff_diff_data);

    upwindSolution->chart()->removeAllSeries();
    laxSolution->chart()->removeAllSeries();
    laxWendroffSolution->chart()->removeAllSeries();
}

void Form::initiateState()
{
    delete param;
    param = new Parameters(spinBoxNX->value()+1, spinBoxNT->value(), kRangeX, kRangeT);

    InitialProfile profile = comboBoxInitial->currentData().value<InitialProfile>();
    state_.resize(param->get_nx()/* + 1*/);
    tmp_state_.resize(state_.size());
    for (decltype(state_.size()) i = 0; i < state_.size(); ++i)
        state_[i] = initial(i * param->get_dx(), profile);

    QList<QPointF> init_data;
    for (decltype(state_.size()) i = 0; i < state_.size(); ++i)
        init_data.append(QPointF(i * param->get_dx(), state_[i]));
    seriesInitial->clear();
    seriesInitial->append(init_data);

    updateLabels();
    updateDispersionDiffusionSolution();
}

void Form::Solve()
{
    pushButtonSolve->setEnabled(false);
    tabWidgetMethods->setEnabled(false);
    comboBoxInitial->setEnabled(false);
    spinBoxNX->setEnabled(false);
    spinBoxNT->setEnabled(false);
    sliderNX->setEnabled(false);
    sliderNT->setEnabled(false);

    initiateState();

    showState();

    t_cur_ = 0.0;
    timer->start();
}

void Form::Tick()
{
    static int t_index = 1;

    if (t_cur_ < kRangeT + 1e-3*param->get_dt())
    {
        t_cur_ += param->get_dt();
        tmp_state_.front() = state_.front();
        tmp_state_.back() = state_.back();
        switch (method_)
        {
        case 0:
            for (decltype(state_.size()) i = 1; i < state_.size()-1; ++i)
                tmp_state_[i] = state_[i] - param->get_alpha() * (state_[i] - state_[i-1]);
            break;
        case 1:
            for (decltype(state_.size()) i = 1; i < state_.size()-1; ++i)
                tmp_state_[i] = 0.5*(state_[i+1] + state_[i-1]) - 0.5*param->get_alpha() * (state_[i+1] - state_[i-1]);
            break;
        case 2:
            for (decltype(state_.size()) i = 1; i < state_.size()-1; ++i)
                tmp_state_[i] = (1.0 - param->get_alpha()*param->get_alpha()) * state_[i] - 0.5*param->get_alpha() * (state_[i+1] - state_[i-1]) + 0.5*param->get_alpha()*param->get_alpha() * (state_[i+1] + state_[i-1]);
            break;
        }

        state_ = tmp_state_;

        if (t_cur_ > kRangeT / 5.0 * t_index)
        {
            ++t_index;
            showState();
        }
    }
    else
    {
        t_index = 1;

        timer->stop();
        pushButtonSolve->setEnabled(true);
        tabWidgetMethods->setEnabled(true);
        comboBoxInitial->setEnabled(true);
        spinBoxNX->setEnabled(true);
        spinBoxNT->setEnabled(true);
        sliderNX->setEnabled(true);
        sliderNT->setEnabled(true);
    }
}

void Form::showState()
{
    QChart *chart = nullptr;
    switch(method_)
    {
    case Upwind:
        chart = upwindSolution->chart();
        break;
    case Lax:
        chart = laxSolution->chart();
        break;
    case LaxWendroff:
        chart = laxWendroffSolution->chart();
        break;
    }

    QLineSeries *series = new QLineSeries();
    chart->addSeries(series);
    series->attachAxis(chart->axisX());
    series->attachAxis(chart->axisY());

    QList<QPointF> data;
    for (decltype(state_.size()) i = 0; i < state_.size(); ++i)
        data << QPointF(i*param->get_dx(), state_[i]);
    series->append(data);
}
