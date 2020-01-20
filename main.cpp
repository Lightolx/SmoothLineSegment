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

// 表示一把弯刀，圆弧再带下面的刀柄
struct SICKLE {
    Eigen::Vector3d q1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d q2 = Eigen::Vector3d::Zero();
    Eigen::Vector3d q3 = Eigen::Vector3d::Zero();
    Eigen::Vector3d O = Eigen::Vector3d::Zero();
    double R = 0.0;
    Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();

    SICKLE() {}
    SICKLE(Eigen::Vector3d _q1,
           Eigen::Vector3d _q2,
           Eigen::Vector3d _q3,
           Eigen::Vector3d _O,
           double _R): q1(std::move(_q1)), q2(std::move(_q2)), q3(std::move(_q3)), O(std::move(_O)), R(_R) {
        ComputeNormal();
    }

    void ComputeNormal() {
        Eigen::Vector3d q2q1 = q1 - q2;
        q2q1.normalize();
        Eigen::Vector3d q2q3 = q3 - q2;
        q2q3.normalize();

        normal = q2q1.cross(q2q3);
        normal.normalize();
    }

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

    // Step1: 生成所有的圆弧
    std::vector<ARC> vArcs;
    for (int i = 1; i < numAnchors - 1; i++) {
//        LOG(INFO) << "i = " << i;
        Eigen::Vector3d p = vAnchors[i];
        Eigen::Vector3d p1 = vAnchors[i - 1];
        Eigen::Vector3d p2 = vAnchors[i + 1];
        
        Eigen::Vector3d t1 = (p1 - p).normalized();
        Eigen::Vector3d t2 = (p2 - p).normalized();

        // 这一段的夹角很小，当成圆弧处理反而误差过大，因为对于几乎是直线的平面算法向量很不稳定
        if (t1.cross(t2).norm() < std::sin(15.0 / 180.0 * M_PI)) {
            continue;
        }

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

        vArcs.emplace_back(q1, q2, O, R);
    }

    // 圆弧再加它后面的那一段直线，共同拼成一把弯刀
    int numArcs = vArcs.size();
    int numSickles = numArcs - 1;
    std::vector<SICKLE> vSickles(numSickles);
    for (int k = 0; k < numArcs-1; ++k) {
        ARC arc1 = vArcs[k];
        ARC arc2 = vArcs[k + 1];
        Eigen::Vector3d q1 = arc1.q1;
        Eigen::Vector3d q2 = arc1.q2;
        Eigen::Vector3d q3 = arc2.q1;

        vSickles[k] = SICKLE(q1, q2, q3, arc1.O, arc1.R);
    }

    std::ofstream fout("../data/Heli_position.txt");
    std::ofstream fout1("../data/Heli_pose.txt");
    for (int k = 0; k < numSickles-1; ++k) {
        LOG(INFO) << "k = " << k;
        SICKLE sickle = vSickles[k];
        Eigen::Vector3d q1 = sickle.q1;
        Eigen::Vector3d q2 = sickle.q2;
        Eigen::Vector3d q3 = sickle.q3;
        Eigen::Vector3d normal = sickle.normal;
        LOG(INFO) << "q1 is " << q1.transpose();
        LOG(INFO) << "q2 is " << q2.transpose();
        LOG(INFO) << "q3 is " << q3.transpose();
        LOG(INFO) << "normal is " << normal.transpose();

        // 对于q1q2, 在刀刃上
        {
            Eigen::Vector3d O = sickle.O;
            double R = sickle.R;
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
                // 计算position
                double alpha = double(i) / R;  // 从r1往r2的方向转，转了alpha的角度
                Eigen::Vector3d n = r1.cross(r2);  // 沿着n轴转动
                n.normalize();
                Eigen::AngleAxisd Ax(alpha, n);
                Eigen::Vector3d r3 = Ax * r1;  // 把r1转alpha角度
                Eigen::Vector3d Oc = O + R * r3;

                // 计算pose
                // 相机的x轴朝右，与r3的方向可能相同也可能相反(r3是从圆心指向车体的向量)，取决于左转还是右转
                // 这里假设car基本上在平面上运动，并且重力方向垂直向下
                Eigen::Vector3d x_axis = Eigen::Vector3d::Zero();
                Eigen::Vector3d y_axis = Eigen::Vector3d::Zero();
                if (normal.dot(Eigen::Vector3d::UnitZ()) > 0) {  // 右转
                    x_axis = -r3;
                    y_axis = -normal;  // 相机的y轴向下
                } else {
                    x_axis = r3;
                    y_axis = normal;  // 相机的y轴向下
                }
                Eigen::Vector3d z_axis = x_axis.cross(y_axis);

                {
                    Eigen::Vector3d l = q2 -q1;
                    l.normalize();
                    if (z_axis.dot(l) < 0) {
                        LOG(ERROR) << "Wrong with the blade direction";
                    }
                }

                Eigen::Matrix3d Rwc = Eigen::Matrix3d::Identity();
                Rwc.col(0) = x_axis;
                Rwc.col(1) = y_axis;
                Rwc.col(2) = z_axis;

                Eigen::Matrix4d Twc = Eigen::Matrix4d::Identity();
                Twc.topLeftCorner(3, 3) = Rwc;
                Twc.topRightCorner(3, 1) = Oc;

                for (int m = 0; m < 3; ++m) {
                    fout << Oc[m] << " ";
                }
                fout << std::endl;

//                LOG(INFO) << numPoints++ << ", in arc, pos = " << Oc.transpose();
                for (int m = 0; m < 3; ++m) {
                    for (int n = 0; n < 4; ++n) {
                        fout1 << Twc(m, n) << " ";
                    }
                }
                fout1 << std::endl;

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

            // 在在直线段还要再加一个功能，就是这个弯刀平面的法向量线性过渡到下一个弯刀平面的法向量
            Eigen::Vector3d normal2 = vSickles[k+1].normal;
            // 计算normal到normal2的夹角，保证均匀变化过去
//            if (normal.dot(normal2) < 0) {  // 说明车在蛇形，刚左转马上又右转
//                normal2 = -normal2;
//            }
            // 虽然normal和normal2可能朝向不一致，也就是一个朝上一个朝下（左转后又右转），
            // 但是无论它俩朝向是否一致，也就是夹得是锐角还是钝角，sin()总是相同的，而且我们想求的是
            // 夹得那个锐角的度数（总不能求个钝角出来吧，那这朝向从normal变化到normal2角度变化也太大了）,
            // 因此这里直接叉乘后再求模长，再求arcsin()函数，刚刚好
            double theta = std::sin(normal.cross(normal2).norm());
            LOG(INFO) << "theta = " << theta;
            int steps = int(length / interval);  // 这个直线段总共要插进去这么多个点
            double interval_angle = theta / steps;  // 每插值一个点，相机的y轴应该从normal向normal2变化这么多角度

            for (double i = 0; i < length; i += interval) {  // 精确到厘米级别
                Eigen::Vector3d Oc = q2 + i * l;

                Eigen::Vector3d y_axis = Eigen::Vector3d::Zero();

                if (normal.dot(Eigen::Vector3d::UnitZ()) > 0) {
                    y_axis = -normal;  // 相机的y轴向下
                } else {
                    y_axis = normal;  // 相机的y轴向下
                }
                Eigen::Vector3d z_axis = l.normalized();
                Eigen::Vector3d x_axis = y_axis.cross(z_axis);
                Eigen::Matrix3d Rwc = Eigen::Matrix3d::Identity();
                Rwc.col(0) = x_axis;
                Rwc.col(1) = y_axis;
                Rwc.col(2) = z_axis;
                Eigen::Matrix4d Twc = Eigen::Matrix4d::Identity();
                Twc.topLeftCorner(3, 3) = Rwc;
                Twc.topRightCorner(3, 1) = Oc;

                for (int m = 0; m < 3; ++m) {
                    fout << Oc[m] << " ";
                }
                fout << std::endl;

//                LOG(INFO) << numPoints++ << ", in line, pos = " << Oc.transpose();
                fout1 << std::fixed << std::setprecision(10);
                for (int m = 0; m < 3; ++m) {
                    for (int n = 0; n < 4; ++n) {
                        fout1 << Twc(m, n) << " ";
                    }
                }
                fout1 << std::endl;

                // normal向量慢慢向normal2向量的方向(或者normal2的反方向)靠拢
                Eigen::Vector3d n = Eigen::Vector3d::Zero();
                if (normal.dot(normal2) > 0) {
                    n = normal.cross(normal2);
                } else {
                    n = normal.cross(-normal2);
                }

                n.normalize();
                // 沿着n轴转动interval_angle度
                Eigen::AngleAxisd Ax(interval_angle, n);
                normal = Ax * normal;
            }
        }
    }

    fout.close();
}