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
#include <hypermesh/ndarray.hh>
#include <hypermesh/regular_simplex_mesh.hh>
#include "trackball.h"

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
 
  void load_ni_trz(const std::string &filename);

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

private:
  hypermesh::ndarray<float> rxy, zxy;

private:
  GLuint tex;
  int current_t = 0;
  int DW, DH, DT;
}; 

#endif
