#ifndef ZCOLOR_H
#define ZCOLOR_H

// credit: Xin Zhang

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


///////////////////////////////
ZRGB::ZRGB()
    : red(0.0), green(0.0), blue(0.0)
{}

ZRGB::ZRGB(double r, double g, double b)
    : red(r), green(g), blue(b)
{}

ZRGB::ZRGB(const ZRGB &src)
    : red(src.red), green(src.green), blue(src.blue)
{}

ZRGB& ZRGB::operator =(const ZRGB &src)
                      {
    red = src.red; green = src.green; blue = src.blue;
    return (*this);
}

void ZRGB::validColor()
{
    if (red > 1.0)  red = 1.0;
    if (red < 0.0) red = 0.0;
    if (green > 1.0) green = 1.0;
    if (green < 0.0) green = 0.0;
    if (blue > 1.0) blue = 1.0;
    if (blue < 0.0) blue = 0.0;
}

ZRGBA::ZRGBA()
    : red(0.0), green(0.0), blue(0.0), alpha(1.0)
{}

ZRGBA::ZRGBA(double r, double g, double b, double a)
    : red(r), green(g), blue(b), alpha(a)
{}

ZRGBA::ZRGBA(const ZRGBA& src)
    : red(src.red), green(src.green), blue(src.blue), alpha(src.alpha)
{}

ZRGBA& ZRGBA::operator=(const ZRGBA& src)
                       {
    red = src.red;
    green = src.green;
    blue = src.blue;
    alpha = src.alpha;
    return (*this);
}

void ZRGBA::validColor()
{
    if (red > 1.0)  red = 1.0;
    if (red < 0.0) red = 0.0;
    if (green > 1.0) green = 1.0;
    if (green < 0.0) green = 0.0;
    if (blue > 1.0) blue = 1.0;
    if (blue < 0.0) blue = 0.0;
    if (alpha > 1.0) alpha = 1.0;
    if (alpha < 0.0) alpha = 0.0;
}

ZLAB::ZLAB()
{}

ZLAB::ZLAB(double l, double a, double b)
    : luminance(l), afactor(a), bfactor(b)
{}

ZLAB::ZLAB(const ZLAB &src)
    : luminance(src.luminance), afactor(src.afactor), bfactor(src.bfactor)
{}

ZLAB& ZLAB::operator =(const ZLAB &src)
                      {
    luminance = src.luminance; afactor = src.afactor; bfactor = src.bfactor;
    return (*this);
}

ZXYZ::ZXYZ()
{}

ZXYZ::ZXYZ(double x, double y, double z)
    : xfactor(x), yfactor(y), zfactor(z)
{}

ZXYZ::ZXYZ(const ZXYZ &src)
    : xfactor(src.xfactor), yfactor(src.yfactor), zfactor(src.zfactor)
{}

ZXYZ& ZXYZ::operator =(const ZXYZ &src)
                      {
    xfactor = src.xfactor; yfactor = src.yfactor; zfactor = src.zfactor;
    return (*this);
}

ZHSI::ZHSI()
{}

ZHSI::ZHSI(double h, double s, double i)
    : hue(h), saturation(s), intensity(i)
{}

ZHSI::ZHSI(const ZHSI &src)
    : hue(src.hue), saturation(src.saturation), intensity(src.intensity)
{}

ZHSI& ZHSI::operator =(const ZHSI &src)
                      {
    hue = src.hue; saturation = src.saturation; intensity = src.intensity;
    return (*this);
}

void ZColor::RGB2XYZ(ZRGB* rgb, ZXYZ* xyz)
{
    if (rgb == 0 || xyz == 0) {
        return;
    }

    double rl, gl, bl;

    if (rgb->red < 0.04045) {
        rl = rgb->red / 12.92;
    } else {
        rl = pow((0.055 + rgb->red) / 1.055, 2.4);
    }
    if (rgb->green < 0.04045) {
        gl = rgb->green / 12.92;
    } else {
        gl = pow((0.055 + rgb->green) / 1.055, 2.4);
    }
    if (rgb->blue < 0.04045) {
        bl = rgb->blue / 12.92;
    } else {
        bl = pow((0.055 + rgb->blue) / 1.055, 2.4);
    }

    xyz->xfactor = 0.4124 * rl + 0.3576 * gl + 0.1805 * bl;
    xyz->yfactor = 0.2126 * rl + 0.7152 * gl + 0.0722 * bl;
    xyz->zfactor = 0.0193 * rl + 0.1192 * gl + 0.9505 * bl;
}

void ZColor::XYZ2LAB(ZXYZ *xyz, ZLAB *lab)
{
    if (xyz == 0 || lab == 0) {
        return;
    }

    double fx, fy, fz;
    double t0 = 6.0 * 6.0 * 6.0 / 29.0 / 29.0 / 29.0;
    double a0 = 29.0 * 29.0 / 6.0 / 6.0 / 3.0;
    double b0 = 4.0 / 29.0;

    if (xyz->xfactor > t0) {
        fx = pow(xyz->xfactor, 1.0 / 3.0);
    } else {
        fx = a0 * xyz->xfactor + b0;
    }
    if (xyz->yfactor > t0) {
        fy = pow(xyz->yfactor, 1.0 / 3.0);
    } else {
        fy = a0 * xyz->yfactor + b0;
    }
    if (xyz->zfactor > t0) {
        fz = pow(xyz->zfactor, 1.0 / 3.0);
    } else {
        fz = a0 * xyz->zfactor + b0;
    }

    lab->luminance = 116.0 * fy - 16.0;
    lab->afactor = 500 * (fx - fy);
    lab->bfactor = 200 * (fy - fz);
}

void ZColor::RGB2LAB(ZRGB *rgb, ZLAB *lab)
{
    if (rgb == 0 || lab == 0) {
        return;
    }

    ZXYZ xyz;
    RGB2XYZ(rgb, &xyz);
    XYZ2LAB(&xyz, lab);
}

void ZColor::LAB2XYZ(ZLAB *lab, ZXYZ *xyz)
{
    if (lab == 0 || xyz == 0) {
        return;
    }

    double fy = (lab->luminance + 16) / 116;
    double fx = fy + lab->afactor / 500;
    double fz = fy - lab->bfactor / 200;
    double delta = 6.0 / 29.0;

    if (fy > delta) {
        xyz->yfactor = fy * fy * fy;
    } else {
        xyz->yfactor = (fy - 16.0 / 116.0) * 3 * delta * delta;
    }
    if (fx > delta) {
        xyz->xfactor = fx * fx * fx;
    } else {
        xyz->xfactor = (fx - 16.0 / 116.0) * 3 * delta * delta;
    }
    if (fz > delta) {
        xyz->zfactor = fz * fz * fz;
    } else {
        xyz->zfactor = (fz - 16.0 / 116.0) * 3 * delta * delta;
    }
}

void ZColor::XYZ2RGB(ZXYZ *xyz, ZRGB *rgb)
{
    if (xyz == 0 || rgb == 0) {
        return;
    }

    double Rl = 3.2410 * xyz->xfactor - 1.5374 * xyz->yfactor - 0.4986 * xyz->zfactor;
    double Gl = -0.9692 * xyz->xfactor + 1.8760 * xyz->yfactor + 0.0416 * xyz->zfactor;
    double Bl = 0.0556 * xyz->xfactor - 0.2040 * xyz->yfactor + 1.0570 * xyz->zfactor;

    if (Rl < 0.00304) {
        rgb->red = 12.92 * Rl;
    } else {
        rgb->red = (1 + 0.055) * pow(Rl, 1.0 / 2.4) - 0.055;
    }
    if (Gl < 0.00304) {
        rgb->green = 12.92 * Gl;
    } else {
        rgb->green = (1 + 0.055) * pow(Gl, 1.0 / 2.4) - 0.055;
    }
    if (Bl < 0.00304) {
        rgb->blue = 12.92 * Bl;
    } else {
        rgb->blue = (1 + 0.055) * pow(Bl, 1.0 / 2.4) - 0.055;
    }

    rgb->validColor();

}

void ZColor::LAB2RGB(ZLAB *lab, ZRGB *rgb)
{
    if (lab == 0 || rgb == 0) {
        return;
    }

    ZXYZ xyz;
    LAB2XYZ(lab, &xyz);
    XYZ2RGB(&xyz, rgb);
}

void ZColor::RGB2HSI(ZRGB *rgb, ZHSI *hsi)
{
    if (rgb == 0 || hsi == 0) {
        return;
    }

    double r = rgb->red;
    double g = rgb->green;
    double b = rgb->blue;
    double m = r;
    if (g < m) m = g;
    if (b < m) m = b;

    double tmpValue = 0.5 * ((r - g) + (r - b)) / sqrt((r - g) * (r - g) + (r - b) * (g - b));
    if (tmpValue > 1.0) {
        tmpValue = 1.0;
    } else if (tmpValue < -1.0) {
        tmpValue = -1.0;
    }
    double theta = acos(tmpValue) / (2 * PI);
    hsi->hue = (b <= g) ? theta : (1 - theta);
    hsi->saturation = 1 - 3 * m / (r + g + b);
    hsi->intensity = (r + g + b) / 3;
}

void ZColor::HSI2RGB(ZHSI *hsi, ZRGB *rgb)
{
    if (hsi == 0 || rgb == 0) {
        return;
    }

    double r, g, b;
    double h = hsi->hue;
    double s = hsi->saturation;
    double i = hsi->intensity;

    h = h * 2 * PI;

    if (h >= 0 && h < 2 * PI / 3)
    {
        b = i * (1 - s);
        r = i * (1 + s * cos(h) / cos(PI / 3 - h));
        g = 3 * i - (r + b);
    } else if (h >= 2 * PI / 3 && h < 4 * PI / 3) {
        r = i * (1 - s);
        g = i * (1 + s * cos(h - 2 * PI / 3) / cos(PI - h));
        b = 3 * i - (r + g);
    } else { //if (h >= 4 * PI / 3 && h <= 2 * PI)
        g = i * (1 - s);
        b = i * (1 + s * cos(h - 4 * PI / 3) / cos(5 * PI / 3 - h));
        r = 3 * i - (g + b);
    }

    rgb->red = r;
    rgb->green = g;
    rgb->blue = b;
    rgb->validColor();
}

#endif
