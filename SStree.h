#ifndef SSTREE_H
#define SSTREE_H

#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include "params.h"
#include "Point.h"
#include <queue>

// ---------------------------------------------------------------------------------------------------------------------

inline int dir_iter = 300;
inline int centroid_iter = 100;

class SsNode {
private:
    [[nodiscard]] static NType varianceAlongDirection(std::vector<Point>::iterator beg,
                                                      std::vector<Point>::iterator end, size_t dim) {
        NType var = 0, mean = 0;
        // VAR[X] = 1/X * sum(X - X_mean)^2

        for (auto it = beg; it != end; ++it) {
            mean += (*it)[dim];
        }

        mean /= (float) std::distance(beg, end);

        for (auto it = beg; it != end; ++it) {
            var += ((*it)[dim] - mean) * ((*it)[dim] - mean);
        }

        var *= 1.0f / (float)std::distance(beg, end);
        return var;
    }

    [[nodiscard]] static size_t minVarianceSplit(std::vector<Point>& centroids, size_t dim) {
        size_t min = Settings::m;
        size_t max = centroids.size() - min;

        NType to_minimize = NType::max_value();
        size_t partition = 0;

        for (size_t i = min; i <= max; ++i) {
            NType var1 = varianceAlongDirection(centroids.begin(), centroids.begin() + i, dim);
            NType var2 = varianceAlongDirection(centroids.begin() + i, centroids.end(), dim);

            if ((var1 + var2) < to_minimize) {
                partition = i;
                to_minimize = var1 + var2;
            }
        }
        return partition;
    }

public:

    [[nodiscard]] static Point meanPoint(const std::vector<Point>& points) {
        Point mean = points[0];

        for (int i = 1; i < points.size(); ++i) {
            mean += points[i];
        }

        mean /= (float) points.size();
        return mean;
    }

    [[nodiscard]] static std::pair<Point, Point> centroidsMBB(const std::vector<Point>& points) {
        size_t dims = points[0].dim();
        Point inf(dims);
        Point sup(dims);

        for (int dim = 0; dim < dims; ++dim) {
            NType min = NType::max_value();
            NType max = NType::min_value();

            for (const Point& p: points) {
                if (p[dim] < min) {
                    min = p[dim];
                }
                if (p[dim] > max) {
                    max = p[dim];
                }
            }

            inf[dim] = min;
            sup[dim] = max;
        }

        return std::make_pair(inf, sup);
    }

public:

    explicit SsNode(std::size_t dims): D(dims) {};
    virtual ~SsNode() = default;

    std::size_t D;
    Point centroid;
    NType radius;
    SsNode* parent = nullptr;

    [[nodiscard]] virtual bool isLeaf() const = 0;
    [[nodiscard]] virtual std::pair<NType, Point> bestCentroidAlongDirection(const Point& p1, const Point& p2, int n) = 0;
    [[nodiscard]] virtual std::vector<Point> getEntriesCentroids() const = 0;
    virtual void sortEntriesByCoordinate(size_t coordinateIndex) = 0;
    virtual std::pair<SsNode*, SsNode*> split() = 0;

    virtual void updateBoundingEnvelope() = 0;

    [[nodiscard]] size_t directionOfMaxVariance() const {
        std::vector<Point> centroids = getEntriesCentroids();
        size_t direction = 0;
        NType max_variance = NType::min_value();

        for (int i = 0; i < centroids[0].dim(); ++i) {
            NType var = varianceAlongDirection(centroids.begin(), centroids.end(), i);

            if (var > max_variance) {
                max_variance = var;
                direction = i;
            }
        }

        return direction;
    }

    [[nodiscard]] size_t findSplitIndex() {
        size_t index = directionOfMaxVariance();
        sortEntriesByCoordinate(index);
        auto centroids = getEntriesCentroids();
        return minVarianceSplit(centroids, index);
    }

    virtual std::pair<SsNode*, SsNode*> insert(const Point& point, const std::string& path) = 0;

    [[nodiscard]] bool test(bool isRoot) const;
    void print(size_t indent = 0) const;
    virtual void saveToStream(std::ostream &out) const = 0;
    virtual void loadFromStream(std::istream &in, SsNode* parent) = 0;
    [[nodiscard]] virtual int getEntriesSize() const = 0;

    [[nodiscard]] virtual int count() const = 0;
    [[nodiscard]] virtual std::pair<NType, Point> FrankWolfeAlgorithm(Point& seed, int k) = 0;
};

// ---------------------------------------------------------------------------------------------------------------------

class SsInnerNode : public SsNode {
private:

    [[nodiscard]] static std::pair<Point, Point> hyperSpheresMBB(const std::vector<SsNode*>& nodes) {
        size_t dims = nodes[0]->centroid.dim();
        Point inf(dims);
        Point sup(dims);

        for (int dim = 0; dim < dims; ++dim) {
            NType min = NType::max_value();
            NType max = NType::min_value();

            for (SsNode* node : nodes) {
                NType current_min = node->centroid[dim] - node->radius;
                NType current_max = node->centroid[dim] + node->radius;

                if (current_min < min) {
                    min = current_min;
                }
                if (current_max > max) {
                    max = current_max;
                }
            }

            inf[dim] = min;
            sup[dim] = max;
        }

        return std::make_pair(inf, sup);
    }

    [[nodiscard]] static std::pair<Point, Point> hyperSpheresSecondMBB(std::pair<Point, Point>& mbb,
                                                                       const std::vector<SsNode*>& nodes) {
        size_t dims = nodes[0]->centroid.dim();
        Point inf(dims);
        Point sup(dims);

        for (int dim = 0; dim < dims; ++dim) {
            NType min = NType::max_value();
            NType max = NType::min_value();

            for (SsNode* node : nodes) {
                NType current_min = node->centroid[dim] - node->radius;
                NType current_max = node->centroid[dim] + node->radius;

                if (current_min < min && current_min > mbb.first[dim]) {
                    min = current_min;
                }
                if (current_max < min && current_max > mbb.first[dim]) {
                    min = current_max;
                }

                if (current_max > max && current_max < mbb.second[dim]) {
                    max = current_max;
                }
                if (current_min > max && current_min < mbb.second[dim]) {
                    max = current_min;
                }
            }

            if (min == NType::max_value()) {
                min = mbb.first[dim];
            }

            if (max == NType::min_value()) {
                max = mbb.second[dim];
            }

            inf[dim] = min;
            sup[dim] = max;
        }

        return std::make_pair(inf, sup);
    }

    [[nodiscard]] static std::pair<Point, Point> hyperSpheresIMBB(const std::vector<SsNode*>& nodes) {
        size_t dims = nodes[0]->centroid.dim();
        Point inf(dims);
        Point sup(dims);

        for (int dim = 0; dim < dims; ++dim) {
            NType min = NType::min_value();
            NType max = NType::max_value();

            for (SsNode* node : nodes) {
                NType current_min = node->centroid[dim] - node->radius;
                NType current_max = node->centroid[dim] + node->radius;

                if (current_min > min) {
                    sup[dim] = current_min;
                }
                if (current_max < max) {
                    inf[dim] = current_max;
                }
            }
        }

        return std::make_pair(inf, sup);
    }

    [[nodiscard]] std::vector<Point> getEntriesCentroids() const override {
        std::vector<Point> centroids;
        centroids.reserve(children.size());

        for (SsNode* node: children) {
            centroids.push_back(node->centroid);
        }

        return centroids;
    }

    void sortEntriesByCoordinate(size_t coordinateIndex) override {
        std::sort(children.begin(), children.end(), [&](SsNode* n1, SsNode* n2) {
            return n1->centroid[coordinateIndex] < n2->centroid[coordinateIndex];
        });
    }

public:

    explicit SsInnerNode(std::size_t dims): SsNode(dims) {}

    std::pair<SsNode*, SsNode*> split() override {
        size_t splitIndex = findSplitIndex();

        SsNode* newNode1 = new SsInnerNode(D);
        dynamic_cast<SsInnerNode*>(newNode1)->children = std::vector<SsNode*>(this->children.begin(),
                                                                              children.begin() + splitIndex);

        SsNode* newNode2 = new SsInnerNode(D);
        dynamic_cast<SsInnerNode*>(newNode2)->children = std::vector<SsNode*>(this->children.begin() + splitIndex,
                                                                              this->children.end());

        newNode1->updateBoundingEnvelope();
        newNode2->updateBoundingEnvelope();

        newNode1->parent = this->parent;
        newNode2->parent = this->parent;

        return  std::make_pair(newNode1, newNode2);
    }

    std::vector<SsNode*> children;

    [[nodiscard]] SsNode* findClosestChild(const Point& target) const {
        SsNode* closest = nullptr;
        NType best_distance = NType::max_value();

        for (SsNode* node: children) {
            NType dis = distance(target, node->centroid);
            if (dis < best_distance) {
                best_distance = dis;
                closest = node;
            }
        }

        return closest;
    }

    [[nodiscard]] SsNode* findMinRadiusIncreaseChild(const Point& target) const {
        SsNode* closest = nullptr;
        NType min_radius_increase = NType::max_value();
        NType num_entries = NType::min_value();

        for (SsNode* node: children) {
            NType dis = distance(node->centroid, target);

            if (dis <= node->radius) {
                if (num_entries < getEntriesSize()) {
                    closest = node;
                    num_entries = getEntriesSize();
                }
            }
            else {
                Point difference = target - node->centroid;
                Point direction = difference / difference.norm();

                Point bounding_point = node->centroid + direction * node->radius;

                NType radius_increase = distance(bounding_point, target);

                if (radius_increase < min_radius_increase) {
                    min_radius_increase = radius_increase;
                    closest = node;
                }
            }
        }

        return closest;
    }

    [[nodiscard]] bool isLeaf() const override { return false; }

    std::pair<NType, Point> bestCentroidAlongDirection(const Point& p1, const Point& p2, int n_iter) override {
        NType segment_distance = distance(p1, p2);
        NType step = segment_distance / (float) n_iter;

        Point difference = p2 - p1;
        Point direction = difference / difference.norm();

        Point best_centroid;           // Initially unassigned
        Point current_centroid = p1;  // Initial centroid

        NType current_centroid_radius; // Initially unassigned
        NType best_radius = NType::max_value(); // Initial best

        for (int i = 0; i < n_iter; ++i){
            current_centroid_radius = 0;
            for (SsNode* node: children) {
                NType dis = distance(node->centroid, current_centroid) + node->radius;
                if (dis > current_centroid_radius) {
                    current_centroid_radius = dis;
                }
            }

            if (current_centroid_radius < best_radius) {
                best_radius = current_centroid_radius;
                best_centroid = current_centroid;
            }

            current_centroid += direction * step;
        }

        return std::make_pair(best_radius, best_centroid);
    }

    void updateBoundingEnvelope() override {
        std::vector<Point> centroids = getEntriesCentroids();

        std::pair<Point, Point> mbb = hyperSpheresMBB(children);
        std::pair<Point, Point> mbb2nd = hyperSpheresSecondMBB(mbb, children);
        std::pair<Point, Point> imbb = hyperSpheresIMBB(children);
        std::pair<Point, Point> mbb_centroids = centroidsMBB(centroids);

        Point mean = meanPoint(centroids);

        std::vector<Point> points = {
                mean,                                               // Promedio de todos los centroides
                (mbb_centroids.first + mbb_centroids.second) / 2.0, // Centro del MBB de centroides
                (mbb.first + mbb.second) / 2.0,                     // Centro del MBB de hiperesferas
                (imbb.first + imbb.second) / 2.0                    // Centro del MBB invertido
        };

        Point reference = meanPoint(points);
        points.push_back(reference);

        NType best_radius = NType::max_value();
        Point best_centroid;

        for (int i = 0; i < points.size(); ++i) {
            for (int j = i; j < points.size(); ++j) {
                if (distance(points[i], points[j]) > 0) {
                    std::pair<NType, Point> option = bestCentroidAlongDirection(points[i], points[j], dir_iter);
                    if (option.first < best_radius) {
                        best_radius = option.first;
                        best_centroid = option.second;
                    }
                }
            }
        }

        if (distance(mbb2nd.first, mbb2nd.second) > 0) {
            std::pair<NType, Point> option = bestCentroidAlongDirection(mbb2nd.first, mbb2nd.second, dir_iter);
            if (option.first < best_radius) {
                best_radius = option.first;
                best_centroid = option.second;
            }
        }

        auto x = FrankWolfeAlgorithm(mean, centroid_iter);
        if (x.first < best_radius) {
            this->centroid = x.second;
            this->radius = x.first;
            return;
        }

        this->radius = best_radius;
        this->centroid = best_centroid;
    }

    void updateBoundingEnvelopeSplit() {
        std::vector<Point> centroids = getEntriesCentroids();

        std::pair<Point, Point> mbb = hyperSpheresMBB(children);
        std::pair<Point, Point> mbb2nd = hyperSpheresSecondMBB(mbb, children);
        std::pair<Point, Point> imbb = hyperSpheresIMBB(children);
        std::pair<Point, Point> mbb_centroids = centroidsMBB(centroids);

        Point mean = meanPoint(centroids);

        std::vector<Point> points = {
                mean,                                               // Promedio de todos los centroides
                (mbb_centroids.first + mbb_centroids.second) / 2.0, // Centro del MBB de centroides
                (mbb.first + mbb.second) / 2.0,                     // Centro del MBB de hiperesferas
                (imbb.first + imbb.second) / 2.0,                   // Centro del MBB invertido
        };

        Point reference = meanPoint(points);
        points.push_back(reference);

        NType best_radius = NType::max_value();
        Point best_centroid;

        for (int i = 0; i < points.size(); ++i) {
            for (int j = i; j < points.size(); ++j) {
                if (distance(points[i], points[j]) > 0) {
                    std::pair<NType, Point> option = bestCentroidAlongDirection(points[i], points[j], dir_iter);
                    if (option.first < best_radius) {
                        best_radius = option.first;
                        best_centroid = option.second;
                    }
                }
            }
        }

        if (distance(mbb2nd.first, mbb2nd.second) > 0) {
            std::pair<NType, Point> option = bestCentroidAlongDirection(mbb2nd.first, mbb2nd.second, dir_iter);
            if (option.first < best_radius) {
                best_radius = option.first;
                best_centroid = option.second;
            }
        }

        if (children.size() > 2) {
            Point c1 = children[children.size() - 1]->centroid;
            NType r1 = children[children.size() - 1]->radius;

            Point c2 = children[children.size() - 2]->centroid;
            NType r2 = children[children.size() - 2]->radius;

            Point v = c1 - c2;
            Point dir = v / v.norm();

            Point lim1 = c1 + dir * r1;
            Point lim2 = c2 - dir * r2;

            Point center = (lim1 + lim2) / 2.0;

            std::vector<Point> newEntriesPoints = {
                    c1,
                    c2,
                    (c1 + c2) / 2.0,
                    center
            };

            for (Point &point: newEntriesPoints) {
                if (distance(centroid, point) > 0) {
                    std::pair<NType, Point> option = bestCentroidAlongDirection(centroid, point, dir_iter);
                    if (option.first < best_radius) {
                        best_radius = option.first;
                        best_centroid = option.second;
                    }
                }
            }
        }

        auto x = FrankWolfeAlgorithm(mean, centroid_iter);
        if (x.first < best_radius) {
            this->centroid = x.second;
            this->radius = x.first;
            return;
        }

        this->radius = best_radius;
        this->centroid = best_centroid;
    }

    [[nodiscard]] std::pair<NType, Point> FrankWolfeAlgorithm(Point& seed, int k) override {
        Point c = seed;
        NType r = NType::max_value();
        int i = 1;
//
//        std::vector<Point> points;
//        for (SsNode* node: children) {
//            Point p = node->centroid;
//            Point v = c - p;
//            Point dir = v / v.norm();
//            points.push_back(c + dir * node->radius);
//        }

        while (true) {
            NType current_radius = NType::min_value();
            Point farthest;
            for (SsNode* node: children) {
                NType d = distance(node->centroid, c);
                if (d > current_radius) {
                    current_radius = d;
                    farthest = node->centroid;
                }
            }

            r = current_radius;

            if (i == k) {
                break;
            }

            c = c * float(i) / float(i + 1) + farthest * 1 / float(i + 1);
            ++i;
        }

        for (SsNode* node: children) {
            NType dis = distance(node->centroid, c) + node->radius;
            if (dis > r) {
                r = dis;
            }
        }

        return {r, c};
    }

    std::pair<SsNode*, SsNode*> insert(const Point& point, const std::string& path) override {
        SsNode* closestChild = findMinRadiusIncreaseChild(point);

        std::pair<SsNode*, SsNode*> newChildren = closestChild->insert(point, path);
        SsNode* newChild1 = newChildren.first;
        SsNode* newChild2 = newChildren.second;

        if (newChild1 == nullptr) {
            updateBoundingEnvelope();
            return std::make_pair(nullptr, nullptr);
        } else {
            children.erase(std::remove(children.begin(),children.end(),closestChild),
                           children.end());
            children.emplace_back(newChild1);
            children.emplace_back(newChild2);
            updateBoundingEnvelopeSplit();

            if (children.size() <= Settings::M) {
                return std::make_pair(nullptr, nullptr);
            }
        }

        return split();
    }

    void saveToStream(std::ostream &out) const override;
    void loadFromStream(std::istream &in, SsNode* parent) override;

    [[nodiscard]] int count() const override {
        int c = 0;
        for (auto& child : children) {
            if (child != nullptr) {
                c += child->count();
            }
        }
        return c;
    }

    [[nodiscard]] int getEntriesSize() const override {
        return this->children.size();
    }
};

// ---------------------------------------------------------------------------------------------------------------------

class SsLeaf : public SsNode {
private:

    [[nodiscard]] std::vector<Point> getEntriesCentroids() const override {
        return points;
    }

    void sortEntriesByCoordinate(size_t coordinateIndex) override {
        int i, j;
        bool swapped;
        size_t n = points.size();

        for (i = 0; i < n - 1; i++) {
            swapped = false;
            for (j = 0; j < n - i - 1; j++) {
                if (points[j][coordinateIndex] > points[j + 1][coordinateIndex]) {
                    std::swap(points[j], points[j + 1]);
                    std::swap(paths[j], paths[j + 1]);
                    swapped = true;
                }
            }
            if (!swapped) {
                break;
            }
        }
    }

public:

    std::vector<Point> points;
    std::vector<std::string> paths;

    explicit SsLeaf(size_t dims): SsNode(dims) {};

    std::pair<SsNode*, SsNode*> split() override {
        size_t splitIndex = findSplitIndex();

        auto* newNode1 = new SsLeaf(D);
        newNode1->points.reserve(splitIndex);
        for (size_t i = 0; i < splitIndex; ++i) {
            newNode1->points.push_back(points[i]);
            newNode1->paths.push_back(paths[i]);
        }
        newNode1->parent = parent;

        auto* newNode2 = new SsLeaf(D);
        newNode2->points.reserve(points.size() - splitIndex);
        for (size_t i = splitIndex; i < points.size(); ++i) {
            newNode2->points.push_back(points[i]);
            newNode2->paths.push_back(paths[i]);
        }
        newNode2->parent = parent;

        newNode1->updateBoundingEnvelope();
        newNode2->updateBoundingEnvelope();
        return std::make_pair(newNode1, newNode2);
    }

    [[nodiscard]] bool isLeaf() const override { return true; }

    std::pair<NType, Point> bestCentroidAlongDirection(const Point& p1, const Point& p2, int n_iter) override {
        NType segment_distance = distance(p1, p2);
        NType step = segment_distance / (float) n_iter;

        Point difference = p2 - p1;
        Point direction = difference / difference.norm();

        Point best_centroid;           // Initially unassigned
        Point current_centroid = p1;   // Initial centroid

        NType current_centroid_radius; // Initially unassigned
        NType best_radius = NType::max_value(); // Initial best

        for (int i = 0; i < n_iter; ++i){
            current_centroid_radius = 0;

            for (Point& point: points) {
                NType dis = distance(point, current_centroid);
                if (dis > current_centroid_radius) {
                    current_centroid_radius = dis;
                }
            }

            if (current_centroid_radius < best_radius) {
                best_radius = current_centroid_radius;
                best_centroid = current_centroid;
            }

            current_centroid += direction * step;
        }

        return std::make_pair(best_radius, best_centroid);
    }

    [[nodiscard]] std::pair<NType, Point> FrankWolfeAlgorithm(Point& seed, int k) override {
        Point c = seed;
        NType r = NType::max_value();
        int i = 1;

        while (true) {
            NType current_radius = NType::min_value();
            Point farthest;
            for (Point &p: points) {
                NType d = distance(p, c);
                if (d > current_radius) {
                    current_radius = d;
                    farthest = p;
                }
            }

            r = current_radius;

            if (i == k) {
                break;
            }

            c = c * float(i) / float(i + 1) + farthest * 1 / float(i + 1);
            ++i;
        }

        return {r, c};
    }

    void updateBoundingEnvelope() override {
        if (points.size() == 1) {
            this->radius = 0;
            this->centroid = points[0];
            return;
        }
        else if (points.size() == 2) {
            this->centroid = (points[0] + points[1]) / 2.0;
            this->radius = distance(points[0], centroid);
            return;
        }

        Point mean = meanPoint(points);
        auto frankWolfeAlgorithm = FrankWolfeAlgorithm(mean, centroid_iter);

        this->centroid = frankWolfeAlgorithm.second;
        this->radius = frankWolfeAlgorithm.first;
    }

    std::pair<SsNode*, SsNode*> insert(const Point& point, const std::string& path) override {
        if (std::find(points.begin(), points.end(), point) != points.end()) {
            return std::make_pair(nullptr, nullptr);
        }

        points.push_back(point);
        paths.push_back(path);

        updateBoundingEnvelope();

        if (points.size() <= Settings::M) {
            return std::make_pair(nullptr, nullptr);
        }

        return split();
    }

    void saveToStream(std::ostream &out) const override;
    void loadFromStream(std::istream &in, SsNode* parent) override;

    [[nodiscard]] int count() const override {
        return points.size();
    }

    [[nodiscard]] int getEntriesSize() const override {
        return points.size();
    }
};

// ---------------------------------------------------------------------------------------------------------------------

class SsTree {
private:

    std::function<bool(std::pair<NType, std::pair<Point, std::string>>&, std::pair<NType, std::pair<Point, std::string>>&)> compare = [&](std::pair<NType, std::pair<Point, std::string>>& a, std::pair<NType, std::pair<Point, std::string>>& b) {
        return a.first < b.first;
    };

    SsNode* root;
    std::size_t D;
//    SsNode* search(SsNode* node, const Point& target);
//    SsNode* searchParentLeaf(SsNode* node, const Point& target);
    void _kNNQuery(const Point &target, size_t k, SsNode *node, NType radius, std::priority_queue<std::pair<NType, std::pair<Point, std::string>>, std::vector<std::pair<NType, std::pair<Point, std::string>>>, decltype(compare)> &result) const {
        if (node->isLeaf()) {
            auto* leafNode = dynamic_cast<SsLeaf*>(node);
            for (int i = 0; i < leafNode->points.size(); ++i) {
                Point p = leafNode->points[i];

                NType dist = std::max(NType(0), distance(target, p));
                // Regla 2
                if (radius + std::max(NType(0), dist) < std::max(NType(0), distance(target, leafNode->centroid))) {
                    continue;
                }
                // Regla 4
                if (radius + std::max(NType(0), distance(target, leafNode->centroid)) < std::max(NType(0), distance(p, leafNode->centroid))) {
                    continue;
                }

                if (dist < radius) {
                    result.push({dist, {p, leafNode->paths[i]}});
                    if (result.size() > k) {
                        result.pop();
                    }
                }
            }
        }
        else {
            auto* innerNode = dynamic_cast<SsInnerNode*>(node);
            // Se obtiene la menor distancia hacia un nodo
            for (auto* child : innerNode->children) {
                NType dist = std::max(NType(0), distance(target, child->centroid));

                // Regla 1
                if (radius + child->radius < std::max(NType(0), distance(target, child->centroid))) {
                    continue;
                }

                // Regla 3
                if (radius + std::max(NType(0), distance(target, child->centroid)) < child->radius) {
                    continue;
                }

                if (dist + child->radius < radius) {
                    _kNNQuery(target, k, child, radius, result);
                }
            }
        }
    }

public:
    explicit SsTree(std::size_t dims) : root(nullptr), D(dims) {}
    ~SsTree() {
        delete root;
    }

    void insert(const Point& point, const std::string& path) {
        if (root == nullptr) {
            root = new SsLeaf(D);
            root->parent = nullptr;
            root->centroid = point;
            root->radius = 0;
        }

        std::pair<SsNode*, SsNode*> newChildren = root->insert(point, path);
        SsNode* newChild1 = newChildren.first;
        SsNode* newChild2 = newChildren.second;

        if (newChild1 == nullptr) {
            return;
        }

        root = new SsInnerNode(D);
        root->parent = nullptr;
        newChild1->parent = root;
        newChild2->parent = root;

        dynamic_cast<SsInnerNode*>(root)->children.emplace_back(newChild1);
        dynamic_cast<SsInnerNode*>(root)->children.emplace_back(newChild2);
        (dynamic_cast<SsInnerNode*>(root))->updateBoundingEnvelopeSplit();
    }

    // Se asume que hay más de k puntos en el árbol
    [[nodiscard]] std::vector<std::string> kNNQuery(const Point &center, size_t k) const {
        if (root == nullptr) return std::vector<std::string>();
        std::priority_queue<std::pair<NType, std::pair<Point, std::string>>, std::vector<std::pair<NType, std::pair<Point, std::string>>>, decltype(compare)> result{compare};
        this->_kNNQuery(center, k, this->root, NType::max_value(), result);
        std::vector<std::string> points;
        for (int i = 0; i < k; ++i) {
//            std::cout << result.top().second.first << "\n";
//            std::cout << result.top().first << "\n";
            points.push_back(result.top().second.second);
            result.pop();
        }
        std::vector<std::string> points_inv;
        std::reverse_copy(points.begin(), points.end(), std::back_inserter(points_inv));
        return points_inv;
    }

    void print() const;
    void test() const;

    void saveToFile(const std::string &filename) const;
    void loadFromFile(const std::string &filename);

    int count() {
        if (!root) {
            return -1;
        }
        return root->count();
    }
};


#endif // !SSTREE_H