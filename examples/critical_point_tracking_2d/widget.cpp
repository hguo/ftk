#include <QMouseEvent>
#include <QFileDialog>
#include <QInputDialog>
#include <QDebug>
#include <fstream>
#include <iostream>
#include <queue>
// #include <netcdf.h>
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

CGLWidget::CGLWidget(const hypermesh::ndarray<float> &s, const QGLFormat& fmt, QWidget *parent, QGLWidget *sharedWidget)
  : QGLWidget(fmt, parent, sharedWidget), 
    scalar(s),
    fovy(30.f), znear(0.1f), zfar(10.f), 
    eye(0, 0, 2.5), center(0, 0, 0), up(0, 1, 0)
{
  DW = s.dim(0);
  DH = s.dim(1);
  DT = s.dim(2);

  auto [min, max] = scalar.min_max();
  scalar_min = min; 
  scalar_max = max;
}

CGLWidget::~CGLWidget()
{
}

void CGLWidget::set_trajectories(const std::vector<std::vector<std::vector<float>>>& traj)
{
  trajectories = traj;
  colors.clear();

  for (int i = 0; i < traj.size(); i ++)
    colors.push_back(QColor::fromHslF((float)rand()/RAND_MAX, 0.5, 0.5));
}

void CGLWidget::mousePressEvent(QMouseEvent* e)
{
  trackball.mouse_rotate(e->x(), e->y()); 
}

void CGLWidget::mouseMoveEvent(QMouseEvent* e)
{
  trackball.motion_rotate(e->x(), e->y()); 
  updateGL(); 
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

  update_texture();
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

  glColor3f(0, 0, 0);
  // glScalef(0.5f, 0.5f, 1.f);
  glTranslatef(-0.5, -0.5, -0.5);
  glPointSize(4.0);
  
  glPushMatrix();
#if 0
  glScalef(1.f/DW, 1.f/DH, 1.f/DT);
  glBegin(GL_POINTS);
  for (const auto &p_ : intersections) {
    const auto &p = p_.second;
    // fprintf(stderr, "%f, %f, %f\n", p.x[0], p.x[1], p.x[2]);
    const auto &c = cc_colors[p.label];
    glColor3f(c.redF(), c.greenF(), c.blueF());
    glVertex3f(p.x[0], p.x[1], p.x[2]);
  }
  glEnd();
#endif

  // glScalef(1.f/DW, 1.f/DH, 1.f/DT);
  glLineWidth(4.0);
  for (int i = 0; i < trajectories.size(); i ++) {
    glColor3f(colors[i].redF(), colors[i].greenF(), colors[i].blueF());
    for (int j = 0; j < trajectories[i].size(); j ++) {
      const auto &c = trajectories[i][j];
      glBegin(GL_LINE_STRIP);
      for (int k = 0; k < c.size()/3; k ++)
        glVertex3f(c[k*3], c[k*3+1], c[k*3+2]);
      glEnd();
    }
  }
  glPopMatrix();

  CHECK_GLERROR();
}

void CGLWidget::update_texture()
{
  hypermesh::ndarray<unsigned char> colors({3, (size_t)DW, (size_t)DH});
  // fprintf(stderr, "current_t=%d\n", current_t);

  for (int j = 0; j < DH; j ++) {
    for (int i = 0; i < DW; i ++) {
      const auto s = (scalar(i, j, current_t) - scalar_min) / (scalar_max - scalar_min) * 255;
      colors(0, i, j) = s; 
      colors(1, i, j) = 255 - s;
      colors(2, i, j) = 0;
    }
  }

  glBindTexture(GL_TEXTURE_2D, tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, DW, DH, 0, GL_RGB, GL_UNSIGNED_BYTE, colors.data());
}

