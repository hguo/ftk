#include "widget.h"
#include <QApplication>
#include <iostream>

int main(int argc, char **argv)
{
  if (argc < 2) return 0;

  QApplication app(argc, argv);
  QGLFormat fmt = QGLFormat::defaultFormat();
  fmt.setSampleBuffers(true);
  fmt.setSamples(16);
  QGLFormat::setDefaultFormat(fmt);

  CGLWidget *widget = new CGLWidget;
  widget->load_ni_trz(argv[1]);
  widget->show();
  return app.exec();
}
