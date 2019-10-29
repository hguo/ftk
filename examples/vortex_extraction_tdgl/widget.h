#ifndef _WIDGET_H
#define _WIDGET_H

// #include <GL/glew.h>
#include <QGLWidget>
#include <QList>
#include <QVector>
#include <QVector3D>
#include <QMatrix4x4>
#include <mutex>
#include <cmath>
#include <set>
#include <ftk/ndarray.hh>
#include <ftk/hypermesh/regular_simplex_mesh.hh>
#include "trackball/trackball.h"
#include "puncture.h"

class QMouseEvent;
class QKeyEvent; 
class QWheelEvent; 

/* 
 * \class   CGLWidget
 * \author  Hanqi Guo
 * \brief   A light-weight Qt-based 3D viewer
*/
class CGLWidget : public QGLWidget
{
  Q_OBJECT

public:
  CGLWidget(const QGLFormat& fmt=QGLFormat::defaultFormat(), QWidget *parent=NULL, QGLWidget *sharedWidget=NULL); 
  ~CGLWidget(); 
  
protected:
  void initializeGL(); 
  void resizeGL(int w, int h); 
  void paintGL();

  void mousePressEvent(QMouseEvent*); 
  void mouseMoveEvent(QMouseEvent*);
  void keyPressEvent(QKeyEvent*); 
  void wheelEvent(QWheelEvent*); 
  
private:
  void update_texture();

private: // camera
  CGLTrackball trackball;
  QMatrix4x4 projmatrix, mvmatrix; 

  const float fovy, znear, zfar; 
  const QVector3D eye, center, up;

public:
  std::vector<QColor> colors;
}; 

#endif
