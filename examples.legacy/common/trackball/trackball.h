#ifndef _GLTRACKBALL_H
#define _GLTRACKBALL_H

#ifdef __APPLE__
 #include <OpenGL/gl.h>
#else
 #include <GL/gl.h>
#endif
#include <QQuaternion>

/*! \brief The trackball class
 *
 * implementation of mouse rotation and translation
 *
 * \sa http://www.opengl.org/wiki/Trackball
 */

class CGLTrackball
{
public:
  CGLTrackball();
  ~CGLTrackball();
	CGLTrackball(const CGLTrackball &trackball);
	CGLTrackball &operator=(const CGLTrackball &trackball);

	void clone(const CGLTrackball &trackball);

  void init();

  void reshape(int w, int h);
  void mouse_rotate(int u, int v);
  void mouse_translate(int u, int v);
  void motion_rotate(int u, int v);
  void motion_translate(int u, int v);
  void joystick(unsigned int mask, int x, int y, int z);
  void wheel(int delta);

  void applyTransform();
  void applyInverseTransform();
  void applyInverseRotation();    // rotation only. 

  void setScale(float factor) {m_scale = factor;} 
  float getScale() const {return m_scale;}
  void  scale(float factor);
  QVector3D getNormDir();
  QQuaternion getRotation() const;
  QVector3D   getTranslation() const;
  QMatrix4x4  getTransform();
	int getViewportWidth() const { return m_width; }
	int getViewportHeight() const { return m_height;}

  void reset();
  void loadStatus(const char *filename);
  void saveStatus(const char *filename);

  const std::string serialize() const; 
  bool unserialize(const std::string&); 

  void setRotation(const QQuaternion &Q);
  void setTranslation(const QVector3D &T); 

  void rotate(const QQuaternion &Q);

  bool alignToAxes(float err=1e-3); 

protected:
  QQuaternion R; 		//!< Rotation
  QQuaternion incR; 	//!< Increase Rotation
  QVector3D T; 		//!< Translation
  GLfloat M[16];

  void m_Quaternion2Matrix();
  static const QQuaternion axisQuaternions[24]; 

private:
  int origin_x, origin_y; 
  int m_x, m_y;
  int m_width, m_height;
  float m_r; 		    //!< min(m_width, m_height)
  float m_scale;
};

#endif
