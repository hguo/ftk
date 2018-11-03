#include <gtest/gtest.h>

#include <ftk/numerics/invmat.hh>

class MatrixTest : public testing::Test {
 public:
  const float kEps = 1e-5;
};

TEST_F(MatrixTest, invmat2) {
  float m[] = {0.72758122, 0.31241218,
               0.84617905, 0.82793148};
  float inv[4];

  float det = ftk::invmat2(m, inv);

  EXPECT_NEAR(0.3380307546139767, det, kEps);

  EXPECT_NEAR( 2.44927856, inv[0], kEps);
  EXPECT_NEAR(-0.9242123,  inv[1], kEps);
  EXPECT_NEAR(-2.50326054, inv[2], kEps);
  EXPECT_NEAR( 2.15241131, inv[3], kEps);
}

TEST_F(MatrixTest, invmat3) {
  float m[] = {0.9102779,  0.44108077, 0.72642273,
               0.39278198, 0.95680469, 0.02683596,
               0.05335823, 0.86960914, 0.43971526};
  float inv[9];

  float det = ftk::invmat3(m, inv);

  EXPECT_NEAR(0.49721770782016816, det, kEps);

  EXPECT_NEAR( 0.79921696, inv[0], kEps);
  EXPECT_NEAR( 0.8804069,  inv[1], kEps);
  EXPECT_NEAR(-1.37406178, inv[2], kEps);
  EXPECT_NEAR(-0.3444775,  inv[3], kEps);
  EXPECT_NEAR( 0.72705064, inv[4], kEps);
  EXPECT_NEAR( 0.52471497, inv[5], kEps);
  EXPECT_NEAR( 0.58427805, inv[6], kEps);
  EXPECT_NEAR(-1.54469698, inv[7], kEps);
  EXPECT_NEAR( 1.40322755, inv[8], kEps);
}

TEST_F(MatrixTest, invmat4) {
  float m[] = {0.37116367, 0.16844887, 0.99227088, 0.71275718,
               0.88786179, 0.19169413, 0.33589513, 0.40073562,
               0.47178621, 0.34809309, 0.46421167, 0.64434074,
               0.9928173,  0.22639279, 0.8466231,  0.48266757};
  float inv[16];

  float det = ftk::invmat4(m, inv);

  EXPECT_NEAR(-0.036813207788038683, det, kEps);

  EXPECT_NEAR(-0.18813969, inv[0],  kEps);
  EXPECT_NEAR( 1.603472,   inv[1],  kEps);
  EXPECT_NEAR(-0.81077231, inv[2],  kEps);
  EXPECT_NEAR( 0.02888778, inv[3],  kEps);
  EXPECT_NEAR(-4.66259281, inv[4],  kEps);
  EXPECT_NEAR(-6.67423021, inv[5],  kEps);
  EXPECT_NEAR( 5.48385738, inv[6],  kEps);
  EXPECT_NEAR( 5.1058445,  inv[7],  kEps);
  EXPECT_NEAR(-0.07998939, inv[8],  kEps);
  EXPECT_NEAR(-2.51478409, inv[9],  kEps);
  EXPECT_NEAR(-0.08469054, inv[10], kEps);
  EXPECT_NEAR( 2.31908296, inv[11], kEps);
  EXPECT_NEAR( 2.71426273, inv[12], kEps);
  EXPECT_NEAR( 4.24332871, inv[13], kEps);
  EXPECT_NEAR(-0.75591577, inv[14], kEps);
  EXPECT_NEAR(-4.45025938, inv[15], kEps);
}
