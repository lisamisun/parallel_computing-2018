#include "googletest/googletest/include/gtest/gtest.h"
#include "lin_algebra.hpp"
#include "generator.hpp"
#include <vector>

double sumVectorElements(std::vector<double> x) {
    double sum = 0;
    for (int i = 0; i < x.size(); i++) {
        sum += x[i];
    }
    return sum;
}

TEST(dotProductTest, dotProduct) {
    std::vector<double> x, y;

    for (int i = 0; i < 5; i++) {
        x.push_back(sin(i));
        y.push_back(cos(i));
    }

    double dotStandart = sin(0) * cos(0) + sin(1) * cos(1) + sin(2) * cos(2) + sin(3) * cos(3) + sin(4) * cos(4);
    double dotStandartromPython = 0.4312188399711047;

    EXPECT_NEAR(dotProduct(x, y), dotStandart, 0.0001);
    EXPECT_NEAR(dotProduct(x, y), dotStandartromPython, 0.0001);
}

TEST(dotProductPar, dotProduct) {
    std::vector<double> x, y;

    for (int i = 0; i < 5; i++) {
        x.push_back(sin(i));
        y.push_back(cos(i));
    }

    ASSERT_EQ(dotProduct_par(x, y), dotProduct(x, y));
}

TEST(linearCombinationTest, linearCombination) {
    std::vector<double> x, y, x_par;
    double a, b;

    for (int i = 0; i < 5; i++) {
        x.push_back(sin(i));
        y.push_back(cos(i));
    }
    a = 0.5;
    b = -2.8;

    double sumStandart = (0.5 * sin(0) + (-2.8) * cos(0)) + \
                         (0.5 * sin(1) + (-2.8) * cos(1)) + \
                         (0.5 * sin(2) + (-2.8) * cos(2)) + \
                         (0.5 * sin(3) + (-2.8) * cos(3)) + \
                         (0.5 * sin(4) + (-2.8) * cos(4));
    double sumStandartFromPython = 2.0220887769933262;

    double normStandart = sqrt((0.5 * sin(0) + (-2.8) * cos(0)) * (0.5 * sin(0) + (-2.8) * cos(0)) + \
                               (0.5 * sin(1) + (-2.8) * cos(1)) * (0.5 * sin(1) + (-2.8) * cos(1)) + \
                               (0.5 * sin(2) + (-2.8) * cos(2)) * (0.5 * sin(2) + (-2.8) * cos(2)) + \
                               (0.5 * sin(3) + (-2.8) * cos(3)) * (0.5 * sin(3) + (-2.8) * cos(3)) + \
                               (0.5 * sin(4) + (-2.8) * cos(4)) * (0.5 * sin(4) + (-2.8) * cos(4)));
    double normStandartFromPython = 4.673799960603955;

    linearCombination(x, y, a, b);

    EXPECT_NEAR(sumVectorElements(x), sumStandart, 0.0001);
    EXPECT_NEAR(sqrt(dotProduct(x, x)), normStandart, 0.0001);

    EXPECT_NEAR(sumVectorElements(x), sumStandartFromPython, 0.0001);
    EXPECT_NEAR(sqrt(dotProduct(x, x)), normStandartFromPython, 0.0001);
}

TEST(linearCombinationPar, linearCombination) {
    std::vector<double> x, y, x_par;
    double a, b;

    for (int i = 0; i < 5; i++) {
        x.push_back(sin(i));
        y.push_back(cos(i));
        x_par.push_back(sin(i));
    }
    a = 0.5;
    b = -2.8;

    linearCombination(x, y, a, b);
    linearCombination_par(x_par, y, a, b);

    ASSERT_EQ(sumVectorElements(x_par), sumVectorElements(x));
    ASSERT_EQ(sqrt(dotProduct(x_par, x_par)), sqrt(dotProduct(x, x)));
}

TEST(matrixVectorProductTest, matrixVectorProduct) {
    RegularGridGenerator generatorA(2, 2, 2);
    CompressedSparseRowMatrix * a = generatorA.generateDoubleCSRMatrix();

    std::vector<double> x;
    for (int i = 0; i < a->size(); i++) {
        x.push_back(sin(i));
    }

    std::vector<double> y(a->size(), 0);

    a->matrixVectorProduct(x, y);

    double sumStandart = 1.1 * (fabs(sin(0+1+1)) + fabs(sin(0+2+1)) + fabs(sin(0+4+1))) * sin(0) + sin(0+1+1) * sin(1) + sin(0+2+1) * sin(2) + sin(0+4+1) * sin(4) + \
                         sin(1+0+1) * sin(0) + 1.1 * (fabs(sin(1+0+1)) + fabs(sin(1+3+1)) + fabs(sin(1+5+1))) * sin(1) + sin(1+3+1) * sin(3) + sin(1+5+1) * sin(5) + \
                         sin(0+2+1) * sin(0) + 1.1 * (fabs(sin(0+2+1)) + fabs(sin(2+3+1)) + fabs(sin(2+6+1))) * sin(2) + sin(2+3+1) * sin(3) + sin(2+6+1) * sin(6) + \
                         sin(3+1+1) * sin(1) + sin(3+2+1) * sin(2) + 1.1 * (fabs(sin(3+1+1)) + fabs(sin(3+2+1)) + fabs(sin(3+7+1))) * sin(3) + sin(3+7+1) * sin(7) + \
                         sin(4+0+1) * sin(0) + 1.1 * (fabs(sin(4+0+1)) + abs(sin(4+5+1)) + fabs(sin(4+6+1))) * sin(4) + sin(4+5+1) * sin(5) + sin(4+6+1) * sin(6) + \
                         sin(5+1+1) * sin(1) + sin(5+4+1) * sin(4) + 1.1 * (fabs(sin(5+1+1)) + fabs(sin(5+4+1)) + fabs(sin(5+7+1))) * sin(5) + sin(5+7+1) * sin(7) + \
                         sin(6+2+1) * sin(2) + sin(6+4+1) * sin(4) + 1.1 * (fabs(sin(6+2+1)) + fabs(sin(6+4+1)) + fabs(sin(6+7+1))) * sin(6) + sin(6+7+1) * sin(7) + \
                         sin(3+7+1) * sin(3) + sin(5+7+1) * sin(5) + sin(6+7+1) * sin(6) + 1.1 * (fabs(sin(3+7+1)) + fabs(sin(5+7+1)) + abs(sin(6+7+1))) * sin(7);
    double sumStandartFromPython =  2.7122538229117903;  

    double normStandart = sqrt((1.1 * (fabs(sin(0+1+1)) + fabs(sin(0+2+1)) + fabs(sin(0+4+1))) * sin(0) + sin(0+1+1) * sin(1) + sin(0+2+1) * sin(2) + sin(0+4+1) * sin(4)) * \
                          (1.1 * (fabs(sin(0+1+1)) + fabs(sin(0+2+1)) + fabs(sin(0+4+1))) * sin(0) + sin(0+1+1) * sin(1) + sin(0+2+1) * sin(2) + sin(0+4+1) * sin(4)) + \
                          (sin(1+0+1) * sin(0) + 1.1 * (fabs(sin(1+0+1)) + fabs(sin(1+3+1)) + fabs(sin(1+5+1))) * sin(1) + sin(1+3+1) * sin(3) + sin(1+5+1) * sin(5)) * \
                          (sin(1+0+1) * sin(0) + 1.1 * (fabs(sin(1+0+1)) + fabs(sin(1+3+1)) + fabs(sin(1+5+1))) * sin(1) + sin(1+3+1) * sin(3) + sin(1+5+1) * sin(5)) + \
                          (sin(0+2+1) * sin(0) + 1.1 * (fabs(sin(0+2+1)) + fabs(sin(2+3+1)) + fabs(sin(2+6+1))) * sin(2) + sin(2+3+1) * sin(3) + sin(2+6+1) * sin(6)) * \
                          (sin(0+2+1) * sin(0) + 1.1 * (fabs(sin(0+2+1)) + fabs(sin(2+3+1)) + fabs(sin(2+6+1))) * sin(2) + sin(2+3+1) * sin(3) + sin(2+6+1) * sin(6)) + \
                          (sin(3+1+1) * sin(1) + sin(3+2+1) * sin(2) + 1.1 * (fabs(sin(3+1+1)) + fabs(sin(3+2+1)) + fabs(sin(3+7+1))) * sin(3) + sin(3+7+1) * sin(7)) * \
                          (sin(3+1+1) * sin(1) + sin(3+2+1) * sin(2) + 1.1 * (fabs(sin(3+1+1)) + fabs(sin(3+2+1)) + fabs(sin(3+7+1))) * sin(3) + sin(3+7+1) * sin(7)) + \
                          (sin(4+0+1) * sin(0) + 1.1 * (fabs(sin(4+0+1)) + abs(sin(4+5+1)) + fabs(sin(4+6+1))) * sin(4) + sin(4+5+1) * sin(5) + sin(4+6+1) * sin(6)) * \
                          (sin(4+0+1) * sin(0) + 1.1 * (fabs(sin(4+0+1)) + abs(sin(4+5+1)) + fabs(sin(4+6+1))) * sin(4) + sin(4+5+1) * sin(5) + sin(4+6+1) * sin(6)) + \
                          (sin(5+1+1) * sin(1) + sin(5+4+1) * sin(4) + 1.1 * (fabs(sin(5+1+1)) + fabs(sin(5+4+1)) + fabs(sin(5+7+1))) * sin(5) + sin(5+7+1) * sin(7)) * \
                          (sin(5+1+1) * sin(1) + sin(5+4+1) * sin(4) + 1.1 * (fabs(sin(5+1+1)) + fabs(sin(5+4+1)) + fabs(sin(5+7+1))) * sin(5) + sin(5+7+1) * sin(7)) + \
                          (sin(6+2+1) * sin(2) + sin(6+4+1) * sin(4) + 1.1 * (fabs(sin(6+2+1)) + fabs(sin(6+4+1)) + fabs(sin(6+7+1))) * sin(6) + sin(6+7+1) * sin(7)) * \
                          (sin(6+2+1) * sin(2) + sin(6+4+1) * sin(4) + 1.1 * (fabs(sin(6+2+1)) + fabs(sin(6+4+1)) + fabs(sin(6+7+1)))  * sin(6) + sin(6+7+1) * sin(7)) + \
                          (sin(3+7+1) * sin(3) + sin(5+7+1) * sin(5) + sin(6+7+1) * sin(6) + 1.1 * (fabs(sin(3+7+1)) + fabs(sin(5+7+1)) + abs(sin(6+7+1))) * sin(7)) * \
                          (sin(3+7+1) * sin(3) + sin(5+7+1) * sin(5) + sin(6+7+1) * sin(6) + 1.1 * (fabs(sin(3+7+1)) + fabs(sin(5+7+1)) + abs(sin(6+7+1))) * sin(7)));
    double normStandartFromPython = 3.3519231655012263;

    EXPECT_NEAR(sumVectorElements(y), sumStandart, 0.4);
    EXPECT_NEAR(sqrt(dotProduct(y, y)), normStandart, 0.4);

    EXPECT_NEAR(sumVectorElements(y), sumStandartFromPython, 0.0001);
    EXPECT_NEAR(sqrt(dotProduct(y, y)), normStandartFromPython, 0.0001);
}

TEST(matrixVectorProductPar, matrixVectorProduct) {
    RegularGridGenerator generatorA(2, 2, 2);
    CompressedSparseRowMatrix * a = generatorA.generateDoubleCSRMatrix();

    std::vector<double> x;
    for (int i = 0; i < a->size(); i++) {
        x.push_back(sin(i));
    }

    std::vector<double> y(a->size(), 0);
    std::vector<double> y_par(a->size(), 0);

    a->matrixVectorProduct(x, y);
    a->matrixVectorProduct_par(x, y_par);

    ASSERT_EQ(sumVectorElements(y_par), sumVectorElements(y));
    ASSERT_EQ(sqrt(dotProduct(y_par, y_par)), sqrt(dotProduct(y, y)));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}