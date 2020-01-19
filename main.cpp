// 平滑折线段，使之处处可导

#include <iostream>
#include <eigen3/Eigen/Eigen>
#include <fstream>
#include <iomanip>
#include <utility>
#include <glog/logging.h>

struct ARC{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Eigen::Vector3d q1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d q2 = Eigen::Vector3d::Zero();
    Eigen::Vector3d O = Eigen::Vector3d::Zero();
    double R = 0.0;

    ARC() {}
    ARC(Eigen::Vector3d _q1, Eigen::Vector3d _q2, Eigen::Vector3d _O, double _R): q1(std::move(_q1)), q2(std::move(_q2)), O(std::move(_O)), R(_R){}
};

std::vector<Eigen::Vector3d> LoadAnchorPts(const std::string &filename) {
    std::ifstream fin(filename);
    std::string line;
    char semi = 0;
    Eigen::Vector3d pt = Eigen::Vector3d::Zero();
    std::vector<Eigen::Vector3d> vPoints;
    while (getline(fin, line)) {
        std::stringstream ss(line);
        ss >> pt(0) >> semi >> pt(1) >> semi >> pt(2);

        vPoints.emplace_back(pt);
    }
    fin.close();

    return vPoints;
}


int main() {
    // 输入点集，要求顺序排列
    std::vector<Eigen::Vector3d> vAnchors = LoadAnchorPts("../data/picking_list.txt");
    int numAnchors = vAnchors.size();
    LOG(INFO) << "numAnchors = " << numAnchors;

    // 首先求出最短的折线段
    double minLength = (vAnchors[1] - vAnchors[0]).norm();
    for (int i = 0; i < numAnchors - 1; ++i) {
        double length = (vAnchors[i+1] - vAnchors[i]).norm();
        if (length < minLength) {
            minLength = length;
        }
    }

    // 这种方法不设置倒角的曲率半径，转而设置弦长
    double chord_length = minLength * 0.4;  // 最大值是0.5, 但不可设置，否则在此处不平滑
    // 设置每两个插值点之间的欧氏距离，这里设置为0.1
    double interval = 0.1;

    int numArcs = numAnchors - 2;
    std::vector<ARC> vArcs(numArcs);
    for (int i = 1; i < numAnchors - 1; i++) {
        LOG(INFO) << "i = " << i;
        Eigen::Vector3d p = vAnchors[i];
        Eigen::Vector3d p1 = vAnchors[i - 1];
        Eigen::Vector3d p2 = vAnchors[i + 1];
        
        Eigen::Vector3d t1 = (p1 - p).normalized();
        Eigen::Vector3d t2 = (p2 - p).normalized();
        Eigen::Vector3d normal = t1.cross(t2);  // 平面的法向量
        normal.normalize();
        Eigen::Vector3d q1 = p + chord_length * t1;
        Eigen::Vector3d q2 = p + chord_length * t2;

        // 求折线段p1/p/p2的内切圆
        Eigen::Vector3d n1 = normal.cross(t1);
        // n1与t2的夹角必然要是锐角
        if (n1.dot(t2) < 0) {
            n1 = -n1;
        }

        Eigen::Vector3d n2 = normal.cross(t2);
        // n2与t1的夹角必然要是锐角
        if (n2.dot(t1) < 0) {
            n2 = -n2;
        }

        // 求射线[q1, n1]与射线[q2, n2]的交点
        // 假设内切圆的半径为R,则有 q1 + R*n1 = q2 + R*n2, 交点为圆心
        // 有方程 R*(n1-n2) = q2-q1, 这是一个超限方程，这里偷个懒，直接求它们模长的比
        double R = (q2 - q1).norm() / (n1 - n2).norm();
        Eigen::Vector3d O = q1 + R * n1;  // 圆心坐标

        vArcs[i-1] = ARC(q1, q2, O, R);
    }

    std::ofstream fout("../data/Heli_position.txt");
    for (int k = 0; k < numArcs - 1; ++k) {
        ARC arc1 = vArcs[k];
        ARC arc2 = vArcs[k + 1];
        Eigen::Vector3d q1 = arc1.q1;
        Eigen::Vector3d q2 = arc1.q2;
        Eigen::Vector3d q3 = arc2.q1;

        // 所以总共分为两段, q1q2为圆弧，q2q3为直线．
        // 对于q1q2
        {
            Eigen::Vector3d O = arc1.O;
            double R = arc1.R;
            Eigen::Vector3d r1 = (q1 - O).normalized();
            Eigen::Vector3d r2 = (q2 - O).normalized();

            // 从r1转到r2，总共转了theta角
            double theta = 0.0;
            if (r1.dot(r2) >= 0) {  // 夹的是锐角
                theta = std::asin(r1.cross(r2).norm());
            } else {
                theta = M_PI - std::asin(r1.cross(r2).norm());
            }

            // 弧长
            double chord = R * theta; // 弧长也是10cm级别
            for (double i = 0; i < chord; i += interval) {  // 精确到0.001弧度
                double alpha = double(i) / R;  // 从r1往r2的方向转，转了alpha的角度
                Eigen::Vector3d n = r1.cross(r2);  // 沿着n轴转动
                n.normalize();
                Eigen::AngleAxisd Ax(alpha, n);
                Eigen::Vector3d r3 = Ax * r1;  // 把r1转alpha角度
                Eigen::Vector3d Oc = O + R * r3;

                for (int m = 0; m < 3; ++m) {
                    fout << Oc[m] << " ";
                }
                fout << std::endl;
            }
        }

        // 对于q2q3
        {
            Eigen::Vector3d q2q3 = q3 - q2;
            double length = q2q3.norm();
            Eigen::Vector3d l = q2q3.normalized();

            for (double i = 0; i < length; i += interval) {  // 精确到厘米级别
                Eigen::Vector3d Oc = q2 + i * l;
                for (int m = 0; m < 3; ++m) {
                    fout << Oc[m] << " ";
                }
                fout << std::endl;
            }
        }
    }

    fout.close();
}