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

template <typename T>
static T clamp(T min, T max, T val)
{
  return std::max(min, std::min(max, val));
}

template <typename T>
static T clamp_normalize(T min, T max, T val)
{
  return (clamp(min, max, val) - min) / (max - min);
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

void CGLWidget::set_mesh(const ftk::simplex_2d_mesh<> &m)
{
  f_vertices.clear();

  for (int i=0; i<m.n(2); i++) {
    int i0 = m.conn[i*3], i1 = m.conn[i*3+1], i2 = m.conn[i*3+2];
    f_vertices.push_back(m.coords[i0*2]);
    f_vertices.push_back(m.coords[i0*2+1]);
    f_vertices.push_back(m.coords[i1*2]);
    f_vertices.push_back(m.coords[i1*2+1]);
    f_vertices.push_back(m.coords[i2*2]);
    f_vertices.push_back(m.coords[i2*2+1]);
  }
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

void CGLWidget::renderSinglePlane()
{
  if (f_vertices.size() == 0) return;

  glColor3f(0, 0, 0);

  glTranslatef(-1.7, 0, 0);
   
#if 0
  if (toggle_labels) {
    glPointSize(5.f);
    glBegin(GL_POINTS);
    for (int i=0; i<m.nNodes; i++) {
      int label = labels[i];
      if (label != 0) {
        QColor c = label_colors[label];
        glColor3ub(c.red(), c.green(), c.blue());
        glVertex2f(m.coords[i*2], m.coords[i*2+1]);
      }
    }
    glEnd();
  } else if (toggle_mesh) {
#endif
  if (1) {
    if (toggle_wireframe) {
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // render wireframe
    }
    else {
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);

    glVertexPointer(2, GL_FLOAT, 0, f_vertices.data());
    glColorPointer(3, GL_FLOAT, 0, f_colors.data() );
    glDrawArrays(GL_TRIANGLES, 0, f_vertices.size()/2);

    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
  }

  CHECK_GLERROR();
}
