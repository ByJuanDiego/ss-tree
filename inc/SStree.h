#ifndef SSTREE_H
#define SSTREE_H

#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include "params.h"
#include "Point.h"
#include <queue>
#include <functional>

// ---------------------------------------------------------------------------------------------------------------------

inline int frankWolfeIterations = 75;

inline Point meanPoint(const std::vector<Point> &points);

inline std::pair<NType, Point> FrankWolfeAlgorithm(std::vector<Point> &points, Point &seed, int k);

// ---------------------------------------------------------------------------------------------------------------------

class SsNode {
private:

    [[nodiscard]] static NType varianceAlongDirection(std::vector<Point>::iterator beg,
                                                      std::vector<Point>::iterator end,
                                                      size_t dim);

    [[nodiscard]] static long minVarianceSplit(std::vector<Point> &centroids, size_t dim);

public:

    std::size_t D;
    Point centroid;
    NType radius;
    SsNode *parent = nullptr;

    explicit SsNode(std::size_t dims) : D(dims) {};

    virtual ~SsNode() = default;

    [[nodiscard]] virtual bool isLeaf() const = 0;

    [[nodiscard]] virtual std::vector<Point> getEntriesCentroids() const = 0;

    virtual void sortEntriesByCoordinate(size_t coordinateIndex) = 0;

    virtual std::pair<SsNode *, SsNode *> split() = 0;

    virtual void updateBoundingEnvelope() = 0;

    [[nodiscard]] size_t directionOfMaxVariance() const;

    [[nodiscard]] long findSplitIndex();

    virtual std::pair<SsNode *, SsNode *> insert(const Point &point, const std::string &path) = 0;

    [[nodiscard]] bool test(bool isRoot) const;

    virtual void saveToStream(std::ostream &out) const = 0;

    virtual void loadFromStream(std::istream &in, SsNode *parent) = 0;

    [[nodiscard]] virtual int getEntriesSize() const = 0;

};

// ---------------------------------------------------------------------------------------------------------------------

class SsInnerNode : public SsNode {
private:

    [[nodiscard]] std::vector<Point> getEntriesCentroids() const override;

    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:

    std::vector<SsNode *> children;

    explicit SsInnerNode(std::size_t dims);

    std::pair<SsNode *, SsNode *> split() override;

    [[nodiscard]] SsNode *findMinRadiusIncreaseChild(const Point &target) const;

    [[nodiscard]] bool isLeaf() const override;

    std::vector<Point> buildReduction();

    void updateBoundingEnvelope() override;

    std::pair<SsNode *, SsNode *> insert(const Point &point, const std::string &path) override;

    void saveToStream(std::ostream &out) const override;

    void loadFromStream(std::istream &in, SsNode *parent) override;

    [[nodiscard]] int getEntriesSize() const override;
};

// ---------------------------------------------------------------------------------------------------------------------

class SsLeaf : public SsNode {
private:

    [[nodiscard]] std::vector<Point> getEntriesCentroids() const override;

    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:

    std::vector<Point> points;
    std::vector<std::string> paths;

    explicit SsLeaf(size_t dims);

    std::pair<SsNode *, SsNode *> split() override;

    [[nodiscard]] bool isLeaf() const override;

    void updateBoundingEnvelope() override;

    std::pair<SsNode *, SsNode *> insert(const Point &point, const std::string &path) override;

    void saveToStream(std::ostream &out) const override;

    void loadFromStream(std::istream &in, SsNode *parent) override;

    [[nodiscard]] int getEntriesSize() const override;
};

// ---------------------------------------------------------------------------------------------------------------------

class SsTree {
private:

    SsNode *root;
    std::size_t D;

    std::function<bool(std::pair<NType, std::pair<Point, std::string>> &,
                       std::pair<NType, std::pair<Point, std::string>> &)> compare = [&](
            std::pair<NType, std::pair<Point, std::string>> &a, std::pair<NType, std::pair<Point, std::string>> &b) {
        return a.first < b.first;
    };

    void kNNQuery(const Point &target, size_t k, SsNode *node, NType radius,
                   std::priority_queue<std::pair<NType, std::pair<Point, std::string>>, std::vector<std::pair<NType, std::pair<Point, std::string>>>, decltype(compare)> &result) const;

public:

    explicit SsTree(std::size_t dims);

    ~SsTree();

    void insert(const Point &point, const std::string &path);

    [[nodiscard]] std::vector<std::string> kNNQuery(const Point &center, size_t k) const;

    void print() const;

    void test() const;

    void saveToFile(const std::string &filename) const;

    void loadFromFile(const std::string &filename);
};


#endif // !SSTREE_H