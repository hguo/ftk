#include <QMouseEvent>
#include <QFileDialog>
#include <QInputDialog>
#include <QDebug>
#include <fstream>
#include <iostream>
#include <queue>
#include <netcdf.h>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include "widget.h"

#ifdef __APPLE__
#include <OpenGL/glu.h>
// #include <GLUT/glut.h>
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

CGLWidget::CGLWidget(const QGLFormat& fmt, QWidget *parent, QGLWidget *sharedWidget)
  : QGLWidget(fmt, parent, sharedWidget), 
    fovy(30.f), znear(0.1f), zfar(10.f), 
    eye(0, 0, 2.5), center(0, 0, 0), up(0, 1, 0)
{
}

CGLWidget::~CGLWidget()
{
}

void CGLWidget::load_ni_trz(const std::string &filename)
{
  rxy.from_netcdf(filename, "RXY");
  zxy.from_netcdf(filename, "ZXY");

  fprintf(stderr, "%d, %d\n", rxy.shape(0), rxy.shape(1));
  exit(1);
}

void CGLWidget::mousePressEvent(QMouseEvent* e)
{
  if (e->buttons() == Qt::LeftButton) {
    if (e->modifiers() == Qt::ShiftModifier) {
      trackball.mouse_translate(e->x(), e->y());
      updateGL();
    } else {
      trackball.mouse_rotate(e->x(), e->y()); 
      updateGL();
    }
  }
}

void CGLWidget::mouseMoveEvent(QMouseEvent* e)
{
  if (e->buttons() == Qt::LeftButton) {
    if (e->modifiers() == Qt::ShiftModifier) {
      trackball.motion_translate(e->x(), e->y()); 
      updateGL();
    } else {
      trackball.motion_rotate(e->x(), e->y());
      updateGL();
    }
  }
}

void CGLWidget::keyPressEvent(QKeyEvent* e)
{
  switch (e->key()) {
  case Qt::Key_Right:
    current_t = (current_t + 1) % DT;
    update_texture();
    updateGL();
    break;
  
  case Qt::Key_Left:
    current_t = (current_t - 1 + DT) % DT;
    update_texture();
    updateGL();
    break;

  default: break;
  }
}

void CGLWidget::wheelEvent(QWheelEvent* e)
{
  trackball.wheel(e->delta());
  updateGL(); 
}

void CGLWidget::initializeGL()
{
  // glewInit();
  trackball.init();

  // opengl smooth rendering
  {
    glEnable(GL_MULTISAMPLE);

    GLint bufs, samples; 
    glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs); 
    glGetIntegerv(GL_SAMPLES, &samples); 

    glEnable(GL_LINE_SMOOTH); 
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST); 
    
    glEnable(GL_POLYGON_SMOOTH); 
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST); 
    
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1, 1);

    glEnable(GL_BLEND); 
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }
  
  // initialze light for tubes
  {
    GLfloat ambient[]  = {0.1, 0.1, 0.1}, 
            diffuse[]  = {0.5, 0.5, 0.5}, 
            specular[] = {0.8, 0.8, 0.8}; 
    GLfloat dir[] = {0, 0, -1}; 
    GLfloat pos[] = {1, 1, 4, 1};
    GLfloat shiness = 100; 

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient); 
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse); 
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular); 
    glLightfv(GL_LIGHT0, GL_POSITION, pos); 
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, dir); 
    
    GLfloat light1_position[] = {-4.0, 4.0, 0.0, 1.0};
    GLfloat light1_spot_direction[] = {1.0, -1.0, 0.0};

    glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, light1_spot_direction);

    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE); 
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE); 

    glEnable(GL_NORMALIZE); 
    glEnable(GL_COLOR_MATERIAL); 
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular); 
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shiness); 
  }
  
  glGenTextures(1, &tex);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
}

void CGLWidget::resizeGL(int w, int h)
{
  trackball.reshape(w, h); 
  glViewport(0, 0, w, h);

  CHECK_GLERROR(); 
}

void CGLWidget::paintGL()
{
  glClearColor(1, 1, 1, 1); 
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
  
  projmatrix.setToIdentity(); 
  projmatrix.perspective(fovy, (float)width()/height(), znear, zfar); 
  mvmatrix.setToIdentity();
  mvmatrix.lookAt(eye, center, up);
  mvmatrix.rotate(trackball.getRotation());
  mvmatrix.scale(trackball.getScale());
  mvmatrix.translate(trackball.getTranslation());

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glLoadMatrixf(projmatrix.data()); 
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  glLoadMatrixf(mvmatrix.data()); 

  glEnable(GL_DEPTH_TEST);

  glColor3f(0, 0, 0);
  // glutWireTeapot(1.0);

  glPushMatrix();
  // glScalef(0.5f, 0.5f, 1.f);
  glScalef((float)(DW-1) / (DH-1), 1.0, 1.0);
  glRotatef(-90, 0, 0, 1);
  glTranslatef(-0.5, -0.5, -0.5);
  glTranslatef(0, 0, (float)current_t/(DT-1));
  glBindTexture(GL_TEXTURE_2D, tex);
  glEnable(GL_TEXTURE_2D);
  glBegin(GL_QUADS);
  glVertex2f(0, 0); glTexCoord2f(0, 0);
  glVertex2f(1, 0); glTexCoord2f(1, 0);
  glVertex2f(1, 1); glTexCoord2f(1, 1);
  glVertex2f(0, 1); glTexCoord2f(0, 1);
  glEnd();
  glDisable(GL_TEXTURE_2D);
  glPopMatrix();

  CHECK_GLERROR();
}

void CGLWidget::update_texture()
{
  // glBindTexture(GL_TEXTURE_2D, tex);
  // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, DW, DH, 0, GL_RGB, GL_UNSIGNED_BYTE, colors.data());
}

