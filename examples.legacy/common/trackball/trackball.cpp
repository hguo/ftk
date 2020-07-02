#include "trackball.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <QMatrix4x4>

#ifdef _PROTOBUF
#include "trackball.pb.h"
#endif

using namespace std;

const QQuaternion CGLTrackball::axisQuaternions[] = { // all possible posures 
  QQuaternion::fromAxisAndAngle( 0, 0, 1, 0), 
  QQuaternion::fromAxisAndAngle( 0, 0, 1, 90), 
  QQuaternion::fromAxisAndAngle( 0, 0, 1, 180), 
  QQuaternion::fromAxisAndAngle( 0, 0, 1, 270), 

  QQuaternion::fromAxisAndAngle( 0, 0,-1, 0), 
  QQuaternion::fromAxisAndAngle( 0, 0,-1, 90), 
  QQuaternion::fromAxisAndAngle( 0, 0,-1, 180), 
  QQuaternion::fromAxisAndAngle( 0, 0,-1, 270), 

  QQuaternion::fromAxisAndAngle( 0, 1, 0, 0), 
  QQuaternion::fromAxisAndAngle( 0, 1, 0, 90), 
  QQuaternion::fromAxisAndAngle( 0, 1, 0, 180), 
  QQuaternion::fromAxisAndAngle( 0, 1, 0, 270), 

  QQuaternion::fromAxisAndAngle( 0,-1, 0, 0), 
  QQuaternion::fromAxisAndAngle( 0,-1, 0, 90), 
  QQuaternion::fromAxisAndAngle( 0,-1, 0, 180), 
  QQuaternion::fromAxisAndAngle( 0,-1, 0, 270), 

  QQuaternion::fromAxisAndAngle( 1, 0, 0, 0), 
  QQuaternion::fromAxisAndAngle( 1, 0, 0, 90), 
  QQuaternion::fromAxisAndAngle( 1, 0, 0, 180), 
  QQuaternion::fromAxisAndAngle( 1, 0, 0, 270), 

  QQuaternion::fromAxisAndAngle(-1, 0, 0, 0), 
  QQuaternion::fromAxisAndAngle(-1, 0, 0, 90), 
  QQuaternion::fromAxisAndAngle(-1, 0, 0, 180), 
  QQuaternion::fromAxisAndAngle(-1, 0, 0, 270), 
};

CGLTrackball::CGLTrackball()
{
  init();
}

CGLTrackball::~CGLTrackball()
{
}
CGLTrackball::CGLTrackball(const CGLTrackball &trackball)
{
  clone(trackball);
}
CGLTrackball &CGLTrackball::operator=(const CGLTrackball &trackball)
{
  R = trackball.getRotation();
  T = trackball.getTranslation();
  m_scale = trackball.getScale();
  m_width = trackball.getViewportWidth();
  m_height = trackball.getViewportHeight();
  m_r = (std::min)(m_width, m_height);
  m_x = m_y = 0;
  m_Quaternion2Matrix();

  return *this;
}
void CGLTrackball::clone(const CGLTrackball &trackball)
{
  R = trackball.getRotation();
  T = trackball.getTranslation();
  m_scale = trackball.getScale();
  m_width = trackball.getViewportWidth();
  m_height = trackball.getViewportHeight();
  m_r = (std::min)(m_width, m_height);
  m_x = m_y = 0;
  m_Quaternion2Matrix();
}

void CGLTrackball::init()
{
  R = QQuaternion(1.0, 0.0, 0.0, 0.0); 	// no rotation
  /*R *= QQuaternion::fromAxisAndAngle(1, 0, 0, +90);
  R *= QQuaternion::fromAxisAndAngle( 0.0, 0.0, 1.0, 180);*/
  R.normalize();
  T = QVector3D(0.0, 0.0, 0.0);           // no translation
  m_x = m_y = 0;
  m_r = 0;
  m_width = m_height = 0;
  m_scale = 1.0;

  m_Quaternion2Matrix();
}

void CGLTrackball::reset()
{
  init();
}

QVector3D CGLTrackball::getNormDir()
{
  return R.vector();
}

void CGLTrackball::scale(float factor)
{
  m_scale *= factor;
}

QQuaternion CGLTrackball::getRotation() const
{
  return R;
}

QVector3D CGLTrackball::getTranslation() const
{
  return T;
}

void CGLTrackball::reshape(int w, int h)
{
  m_width = w;
  m_height = h;

  m_r = min(m_width, m_height) / 2.f;
}

void CGLTrackball::mouse_rotate(int u, int v)
{
  v = m_height - v;

  origin_x = u; 
  origin_y = v; 

  m_x = u;
  m_y = v;
}

void CGLTrackball::mouse_translate(int u, int v)
{
  v = m_height - v;

  m_x = u;
  m_y = v;
}

void CGLTrackball::motion_rotate(int u, int v)
{
  QVector3D N;
  QVector3D v0(0, 0, 0), v1(0, 0, 0);
  float theta;

  v = m_height - v;

  if (u<0 || u>m_width || v<0 || v>m_height) return;

#if 0
  v0.setX(m_x - m_width*0.5f);
  v0.setY(m_y - m_height*0.5f);
#else
  v0.setX(m_x - origin_x);
  v0.setY(m_y - origin_y);
#endif
  v0.setZ(0.f);

  if (v0.length() < 1e-2) {m_x=u; m_y=v; return;}
  v0.setZ(m_r*m_r/2.f / v0.length());

#if 0
  v1.setX(u - m_width*0.5f);
  v1.setY(v - m_height*0.5f);
#else 
  v1.setX(u - origin_x);
  v1.setY(v - origin_y);
#endif 
  v1.setZ(0.f);

  if (v1.length() < 2) return;
  v1.setZ(m_r*m_r/2.f / v1.length());

  v0.normalize();
  v1.normalize();
  N = QVector3D::crossProduct(v0, v1);
  N.normalize();

  float dotProduct = QVector3D::dotProduct(v0, v1);

  if(dotProduct > 1.0)
    dotProduct = 1.0;

  if(dotProduct < -1.0)
    dotProduct = -1.0;

  theta = acos(dotProduct);

  incR = QQuaternion(cos(theta/2.f), N*sin(theta/2.f));

  R = incR*R;

  m_Quaternion2Matrix();

  m_x = u, m_y = v;
}

void CGLTrackball::motion_translate(int u, int v)
{
  v = m_height - v;

  // a straight forward scheme
  QVector3D screenMovementVector = QVector3D((float)(u - m_x)/(float)m_width, (float)(v - m_y)/(float)m_height, 0.0);

  m_Quaternion2Matrix();
  //  screenMovementVector *= m_scale;

  QMatrix4x4 inverseMatrix = QMatrix4x4(M);
  // QVector3D movementVector = inverseMatrix * screenMovementVector;
  QVector3D movementVector = screenMovementVector;

  T += movementVector;

  m_x = u, m_y = v;
}

void CGLTrackball::rotate(const QQuaternion &Q)
{
  R = Q*R;

  m_Quaternion2Matrix();
}

void CGLTrackball::setRotation(const QQuaternion &Q)
{
  R = Q;

  m_Quaternion2Matrix();
}

void CGLTrackball::setTranslation(const QVector3D &T_)
{
  T = T_;  
}

void CGLTrackball::wheel(int delta)
{
  if (delta > 0)
    m_scale *= 1.05f;
  else
    m_scale *= 0.95f;
}

void CGLTrackball::applyTransform()
{
  m_Quaternion2Matrix();
  glMultMatrixf(M);
  glScalef(m_scale, m_scale, m_scale);
  glTranslatef(T.x(), T.y(), T.z());
}

QMatrix4x4 CGLTrackball::getTransform()
{
  QMatrix4x4 mat;
  mat.setToIdentity();

  mat.rotate(R);
  mat.scale(m_scale);
  mat.translate(T);

  return mat;
}

void CGLTrackball::applyInverseTransform()
{
  QQuaternion S = R;
  R.setScalar(-R.scalar());
  m_Quaternion2Matrix();
  R = S;

  glMultMatrixf(M);
  glScalef(1.f/m_scale, 1.f/m_scale, 1.f/m_scale);
}

void CGLTrackball::applyInverseRotation()
{
  QQuaternion S = R;
  R.setScalar(-R.scalar());
  m_Quaternion2Matrix();
  R = S;

  glMultMatrixf(M);
}

void CGLTrackball::m_Quaternion2Matrix()
{
#if 1
  QMatrix4x4 mat; 
  mat.setToIdentity();
  mat.rotate(R); 
  memcpy(M, mat.data(), sizeof(GLfloat)*16); 
#else
  float q0 = R.scalar(),
      q1 = R.x(),
      q2 = R.y(),
      q3 = R.z();

  M[0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
  M[1] = 2*(q1*q2 + q0*q3);
  M[2] = 2*(q1*q3 - q0*q2);
  M[3] = 0;

  M[4] = 2*(q1*q2 - q0*q3);
  M[5] = q0*q0 + q2*q2 - q1*q1 - q3*q3;
  M[6] = 2*(q2*q3 + q0*q1);
  M[7] = 0;

  M[8] = 2*(q1*q3 + q0*q2);
  M[9] = 2*(q2*q3 - q0*q1);
  M[10]= q0*q0 + q3*q3 - q1*q1 - q2*q2;
  M[11]= 0;

  M[12]= 0;
  M[13]= 0;
  M[14]= 0;
  M[15]= 1;
#endif
}

bool CGLTrackball::alignToAxes(float err)
{
  R.normalize(); 
  for (int i=0; i<24; i++) 
    if ((R - axisQuaternions[i]).lengthSquared() < err) {
      R = axisQuaternions[i]; 
      m_Quaternion2Matrix(); 
      qDebug() << "aligned with" << i; 
      return true; 
    }
  return false; 
}

void CGLTrackball::loadStatus(const char *filename)
{
  ifstream ifs;
  float x, y, z, w, s;

  ifs.open(filename);

  if (!ifs) {
    cerr << "[CGLTrackball::loadStatus] cannot open trackball file. " << endl;
    return;
  }

  ifs >> x >> y >> z >> w >> s;
  R = QQuaternion(w, x, y, z);
  m_scale = s;

  ifs >> x >> y >> z;
  T = QVector3D(x, y, z);

  ifs.close();

  m_Quaternion2Matrix();
}

void CGLTrackball::saveStatus(const char *filename)
{
  ofstream ofs;

  ofs.open(filename);

  if (!ofs) {
    cerr << "[CGLTrackball::saveStatus] cannot save trackball file. " << endl;
    return;
  }

  ofs << R.x() << '\t' << R.y() << '\t' << R.z() << '\t' << R.scalar() << '\t' << m_scale << '\t'
    << T.x() << '\t' << T.y() << '\t' << T.z();
  ofs.close();
}

#ifdef _PROTOBUF
const string CGLTrackball::serialize() const
{
  string buf; 
  PBTrackball pb; 
  pb.set_f_scale(m_scale); 
  pb.add_q_r(R.x()); 
  pb.add_q_r(R.y()); 
  pb.add_q_r(R.z()); 
  pb.add_q_r(R.scalar()); 
  pb.add_v_trans(T.x()); 
  pb.add_v_trans(T.y()); 
  pb.add_v_trans(T.z()); 
  pb.SerializeToString(&buf); 
  return buf; 
}

bool CGLTrackball::unserialize(const string& buf)
{
  PBTrackball pb;  
  if (!pb.ParseFromString(buf)) return false; 
  if (pb.q_r_size()<4 || pb.v_trans_size()<3) return false; 
  m_scale = pb.f_scale(); 
  R.setX(pb.q_r(0)); 
  R.setY(pb.q_r(1)); 
  R.setZ(pb.q_r(2)); 
  R.setScalar(pb.q_r(3)); 
  T.setX(pb.v_trans(0)); 
  T.setY(pb.v_trans(1)); 
  T.setZ(pb.v_trans(2)); 
  return true; 
}
#endif
