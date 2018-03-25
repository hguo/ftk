#ifndef ZCOLOR_H
#define ZCOLOR_H

//sRGB color space, after gamma correlation
//range 0 ~ 1.0 for red, green, blue
class ZRGB
{
public:
        double red;
        double green;
        double blue;
        ZRGB();
        ZRGB(double r, double g, double b);
        ZRGB(const ZRGB& src);
        ZRGB& operator=(const ZRGB& src);
        void validColor();
};

class ZRGBA
{
public:
        double red;
        double green;
        double blue;
        double alpha;
        ZRGBA();
        ZRGBA(double r, double g, double b, double a);
        ZRGBA(const ZRGBA& src);
        ZRGBA& operator=(const ZRGBA& src);
        void validColor();
};

//CIELAB color space
//range 0 ~ 100 for luminance
//different range for afactor and bfactor, aroung -200 ~ 200
class ZLAB
{
public:
        double luminance;
        double afactor;
        double bfactor;
        ZLAB();
        ZLAB(double l, double a, double b);
        ZLAB(const ZLAB& src);
        ZLAB& operator=(const ZLAB& src);
};

//CIEXYZ color space
//different range for xfactor, yfactor, zfactor, around -200 ~ 200
class ZXYZ
{
public:
        double xfactor;
        double yfactor;
        double zfactor;
        ZXYZ();
        ZXYZ(double x, double y, double z);
        ZXYZ(const ZXYZ& src);
        ZXYZ& operator=(const ZXYZ& src);
};

//hsi color space
//range 0 ~ 1.0 for hue, saturation, intensity
class ZHSI
{
public:
        double hue;
        double saturation;
        double intensity;
        ZHSI();
        ZHSI(double h, double s, double i);
        ZHSI(const ZHSI& src);
        ZHSI& operator=(const ZHSI& src);
};

class ZColor
{
public:
        static void RGB2XYZ (ZRGB* rgb, ZXYZ* xyz);
        static void XYZ2LAB (ZXYZ* xyz, ZLAB* lab);
        static void RGB2LAB (ZRGB* rgb, ZLAB* lab);

        static void LAB2XYZ (ZLAB* lab, ZXYZ* xyz);
        static void XYZ2RGB (ZXYZ* xyz, ZRGB* rgb);
        static void LAB2RGB (ZLAB* lab, ZRGB* rgb);

        static void RGB2HSI (ZRGB* rgb, ZHSI* hsi);
        static void HSI2RGB (ZHSI* hsi, ZRGB* rgb);
};

#endif
