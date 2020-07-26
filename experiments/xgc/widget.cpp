#include <QMouseEvent>
#include <QFileDialog>
#include <QInputDialog>
#include <QDebug>
#include <fstream>
#include <iostream>
#include <queue>
#include <functional>
#include <hdf5.h>
// #include <json.hpp>
#include "widget.h"

#ifdef __APPLE__
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glu.h>
// #include <GL/glut.h>
#endif

#define CHECK_GLERROR()\
{\
  GLenum err = glGetError();\
  if (err != GL_NO_ERROR) {\
    const GLubyte *errString = gluErrorString(err);\
    qDebug("[%s line %d] GL Error: %s\n",\
            __FILE__, __LINE__, errString);\
  }\
}

CGLWidget::CGLWidget(const QGLFormat& fmt, QWidget *parent, QGLWidget *sharedWidget) :
  _fovy(30.f), _znear(0.1f), _zfar(10.f), 
  _eye(0, 0, 2.5), _center(0, 0, 0), _up(0, 1, 0)
{
}

CGLWidget::~CGLWidget()
{
}

#if 0
int CGLWidget::advanceSlice() {
  currentSlice ++;
  if (currentSlice >= nphi) currentSlice -= nphi;
  return currentSlice;
}

int CGLWidget::recedeSlice() {
  currentSlice --;
  if (currentSlice < 0) currentSlice += nphi;
  return currentSlice;
}
#endif

int CGLWidget::advanceTimestep() {
  if (currentTimestep < nt-1) currentTimestep ++;
  return currentTimestep;
}

int CGLWidget::recedeTimestep() {
  if (currentTimestep > 0) currentTimestep --;
  return currentTimestep;
}

void CGLWidget::mousePressEvent(QMouseEvent* e)
{
  _trackball.mouse_rotate(e->x(), e->y()); 
}

void CGLWidget::mouseMoveEvent(QMouseEvent* e)
{
  _trackball.motion_rotate(e->x(), e->y()); 
  updateGL(); 
}

void CGLWidget::keyPressEvent(QKeyEvent* e)
{
#if 0
  switch (e->key()) {
  case Qt::Key_Up: 
    advanceTimestep();
    fprintf(stderr, "%d.%d\n", currentTimestep, currentSlice);
    updateDataGL();
    updateGL();
    break;

  case Qt::Key_Down:
    recedeTimestep();
    fprintf(stderr, "%d.%d\n", currentTimestep, currentSlice);
    updateDataGL();
    updateGL();
    break;

  case Qt::Key_Left:
    recedeSlice();
    fprintf(stderr, "%d.%d\n", currentTimestep, currentSlice);
    updateDataGL();
    updateGL();
    break;

  case Qt::Key_Right:
    advanceSlice();
    fprintf(stderr, "%d.%d\n", currentTimestep, currentSlice);
    updateDataGL();
    updateGL();
    break;

  default: break;
  }
#endif
}

void CGLWidget::wheelEvent(QWheelEvent* e)
{
  _trackball.wheel(e->delta());
  updateGL(); 
}

void CGLWidget::initializeGL()
{
  glewInit();
  CHECK_GLERROR();
}

void CGLWidget::resizeGL(int w, int h)
{
  _trackball.reshape(w, h); 
  glViewport(0, 0, w, h);

  CHECK_GLERROR(); 
}

void CGLWidget::paintGL()
{
  glClearColor(1, 1, 1, 0); 
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

  _projmatrix.setToIdentity(); 
  _projmatrix.perspective(_fovy, (float)width()/height(), _znear, _zfar); 
  _mvmatrix.setToIdentity();
  _mvmatrix.lookAt(_eye, _center, _up);
  _mvmatrix.rotate(_trackball.getRotation());
  _mvmatrix.scale(_trackball.getScale());

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glLoadMatrixf(_projmatrix.data()); 
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  glLoadMatrixf(_mvmatrix.data()); 

  glScalef(0.5004, 1, 1);
  glRotatef(-90, 0, 0, 1);
  glTranslatef(-0.5, -0.5, 0);
  glColor4f(1, 1, 1, 1);

  // render mesh
  {

  }

  CHECK_GLERROR();
}

