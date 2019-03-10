#ifndef _FTK_RANDOM_COLOR_H
#define _FTK_RANDOM_COLOR_H

#include <vector>
#include <string>

namespace ftk {

void generate_random_colors(int count, std::vector<unsigned char>& colors);
void generate_colors(int count, std::vector<unsigned char>& colors);


/////////////////////////
inline void generate_random_colors(int count, std::vector<unsigned char>& colors)
{
  const double saturation = 0.7, 
               intensity = 0.5;

  std::vector<double> hues(count), intensities(count);
  for (int i=0; i<count; i++) {
    hues[i] = (double)i/(count-1);
    intensities[i] = (double)i/(count-1)*0.1+0.4;
  }

  std::random_shuffle(hues.begin(), hues.end());
  std::random_shuffle(hues.begin(), hues.end());
  std::random_shuffle(hues.begin(), hues.end());
  std::random_shuffle(hues.begin(), hues.end());
  std::random_shuffle(hues.begin(), hues.end());
  std::random_shuffle(hues.begin(), hues.end());
  std::random_shuffle(hues.begin(), hues.end());
  std::random_shuffle(hues.begin(), hues.end());
  std::random_shuffle(hues.begin(), hues.end());
  std::random_shuffle(hues.begin(), hues.end());
  std::random_shuffle(hues.begin(), hues.end());
  
  std::random_shuffle(intensities.begin(), intensities.end());

  for (int i=0; i<count; i++) {
    ZHSI hsi;
    hsi.hue = hues[i];
    hsi.saturation = saturation;
    hsi.intensity = intensity;
    // hsi.intensity = intensities[i];

    ZRGB rgb;
    ZColor::HSI2RGB(&hsi, &rgb);

    colors.push_back(rgb.red * 255);
    colors.push_back(rgb.green * 255);
    colors.push_back(rgb.blue * 255);
  }
}

inline void generate_colors(int count, std::vector<unsigned char>& colors)
{
  const double saturation = 0.7, 
               intensity = 0.5;

  std::vector<double> hues(count);
  for (int i=0; i<count; i++) 
    hues[i] = (double)i/(count-1);
  
  for (int i=0; i<count; i++) {
    ZHSI hsi;
    hsi.hue = hues[i];
    hsi.saturation = saturation;
    hsi.intensity = intensity;

    ZRGB rgb;
    ZColor::HSI2RGB(&hsi, &rgb);

    colors.push_back(rgb.red * 255);
    colors.push_back(rgb.green * 255);
    colors.push_back(rgb.blue * 255);
  }
}


}

#endif
