#include "form.h"
#include <QApplication>
#include <QTranslator>


int main(int argc, char *argv[])
{   
    QApplication a(argc, argv);

    QTranslator translator;
    translator.load("TransferEquation1D_rus", a.applicationDirPath());
    a.installTranslator(&translator);

    QFont defaultFont;
    defaultFont.setPixelSize(14);
    a.setFont(defaultFont);

    Form w;
    w.show();

    return a.exec();
}
