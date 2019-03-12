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
  CGLWidget(const hypermesh::ndarray<float> &scalar, const QGLFormat& fmt=QGLFormat::defaultFormat(), QWidget *parent=NULL, QGLWidget *sharedWidget=NULL); 
  ~CGLWidget(); 

  void update_texture();

protected:
  void initializeGL(); 
  void resizeGL(int w, int h); 
  void paintGL();

  void mousePressEvent(QMouseEvent*); 
  void mouseMoveEvent(QMouseEvent*);
  void keyPressEvent(QKeyEvent*); 
  void wheelEvent(QWheelEvent*); 

private: // camera
  CGLTrackball trackball;
  QMatrix4x4 projmatrix, mvmatrix; 

  const float fovy, znear, zfar; 
  const QVector3D eye, center, up;

private:
  GLuint tex;
  int current_t = 0;
  int DW, DH, DT;

  const hypermesh::ndarray<float>& scalar;
  float scalar_min, scalar_max;

#if 0
  hypermesh::regular_simplex_mesh m;
  std::map<hypermesh::regular_simplex_mesh_element, intersection_t> intersections;
  std::mutex intersections_mutex;

  std::vector<QColor> cc_colors;
  std::vector<std::vector<std::vector<float>>> curves;

  std::vector<std::tuple<std::vector<unsigned int>, std::vector<float>, std::vector<float>, std::vector<float>>> tubes;
#endif
}; 

#endif
