#include "SStree.h"


NType SsNode::varianceAlongDirection(std::vector<Point>::iterator beg, std::vector<Point>::iterator end, size_t dim) {
    // VAR[X] = 1/X * sum(X - X_mean)^2
    NType var = 0, mean = 0;

    for (auto it = beg; it != end; ++it) {
        mean += (*it)[dim];
    }

    mean /= (float) std::distance(beg, end);

    for (auto it = beg; it != end; ++it) {
        var += ((*it)[dim] - mean) * ((*it)[dim] - mean);
    }

    var *= 1.0f / (float) std::distance(beg, end);
    return var;
}

long SsNode::minVarianceSplit(std::vector<Point> &centroids, size_t dim) {
    size_t min = Settings::m;
    size_t max = centroids.size() - min;

    NType to_minimize = NType::max_value();
    long partition = 0;

    for (long i = (long) min; i <= max; ++i) {
        NType var1 = varianceAlongDirection(centroids.begin(), std::next(centroids.begin(), i), dim);
        NType var2 = varianceAlongDirection(std::next(centroids.begin(), i), centroids.end(), dim);

        if ((var1 + var2) < to_minimize) {
            partition = i;
            to_minimize = var1 + var2;
        }
    }
    return partition;
}

size_t SsNode::directionOfMaxVariance() const {
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

long SsNode::findSplitIndex() {
    size_t index = directionOfMaxVariance();
    sortEntriesByCoordinate(index);
    auto centroids = getEntriesCentroids();
    return minVarianceSplit(centroids, index);
}

bool SsNode::test(bool isRoot) const {
    size_t count;

    if (this->isLeaf()) {
        const SsLeaf *leaf = dynamic_cast<const SsLeaf *>(this);
        count = leaf->points.size();

        // Verificar si los puntos están dentro del radio del nodo
        for (const Point &point: leaf->points) {
            if (distance(this->centroid, point) > this->radius) {
                std::cout << "Point outside node radius detected." << std::endl;
                return false;
            }
        }
    } else {
        const SsInnerNode *inner = dynamic_cast<const SsInnerNode *>(this);
        count = inner->children.size();

        // Verificar si los centroides de los hijos están dentro del radio del nodo padre
        for (const SsNode *child: inner->children) {
            if (distance(this->centroid, child->centroid) > this->radius) {
                std::cout << "Child centroid outside parent radius detected." << std::endl;
                return false;
            }
            // Verificar recursivamente cada hijo
            if (!child->test(false)) {
                return false;
            }
        }
    }

    // Comprobar la validez de la cantidad de hijos/puntos
    if (!isRoot && (count < Settings::m || count > Settings::M)) {
        std::cout << "Invalid number of children/points detected." << std::endl;
        return false;
    }

    // Comprobar punteros de parentezco, salvo para el nodo raíz
    if (!isRoot && !parent) {
        std::cout << "Node without parent detected." << std::endl;
        return false;
    }

    return true;
}

void SsTree::test() const {
    if (!root) {
        std::cout << "SS-Tree is empty!" << std::endl;
    }

    bool result = root->test(true);

    if (root->parent) {
        std::cout << "Root node parent pointer is not null!" << std::endl;
        result = false;
    }

    if (result) {
        std::cout << "SS-Tree is valid!" << std::endl;
    } else {
        std::cout << "SS-Tree has issues!" << std::endl;
    }
}

void SsTree::print() const {
    if (root) {
        std::cout << "Root: " << root->radius << std::endl;

        if (!root->isLeaf()) {
            SsInnerNode *r = dynamic_cast<SsInnerNode *>(root);
            for (auto *x: r->children) {
                std::cout << "          " << x->radius << std::endl;
            }
        }
    } else {
        std::cout << "Empty tree." << std::endl;
    }
}

void SsLeaf::saveToStream(std::ostream &out) const {
    // Guardar centroid
    centroid.saveToFile(out, D);

    // Guardar el radio
    float radius_ = radius.getValue();
    out.write(reinterpret_cast<const char *>(&radius_), sizeof(radius_));

    // Guardar el numero de puntos
    size_t numPoints = points.size();
    out.write(reinterpret_cast<const char *>(&numPoints), sizeof(numPoints));

    // Guardar los puntos
    for (const auto &point: points) {
        point.saveToFile(out, D);
    }

    // Guardar las rutas (paths)
    size_t numPaths = paths.size();
    out.write(reinterpret_cast<const char *>(&numPaths), sizeof(numPaths));
    for (const auto &p: paths) {
        size_t pathLength = p.size();
        out.write(reinterpret_cast<const char *>(&pathLength), sizeof(pathLength));
        out.write(p.c_str(), (long) pathLength);
    }
}

void SsInnerNode::saveToStream(std::ostream &out) const {
    // Guardar centroid
    centroid.saveToFile(out, D);

    // Guardar el radio
    float radius_ = radius.getValue();
    out.write(reinterpret_cast<const char *>(&radius_), sizeof(radius_));

    // Guardar si apunta a nodos hoja
    bool pointsToLeafs = children[0]->isLeaf();
    out.write(reinterpret_cast<const char *>(&pointsToLeafs), sizeof(pointsToLeafs));

    // Guardar la cantidad de hijos para saber cuántos nodos leer después
    size_t numChildren = children.size();
    out.write(reinterpret_cast<const char *>(&numChildren), sizeof(numChildren));

    // Guardar los hijos
    for (const auto &child: children) {
        child->saveToStream(out);
    }
}

void SsInnerNode::loadFromStream(std::istream &in, SsNode *parent) {
    this->parent = parent;

    // Leer centroid
    centroid.readFromFile(in, D);

    // leer el valor del radio
    float radius_ = 0;
    in.read(reinterpret_cast<char *>(&radius_), sizeof(radius_));
    this->radius = radius_;

    // leer si apunta a hojas o nodos internos
    bool pointsToLeaf = false;
    in.read(reinterpret_cast<char *>(&pointsToLeaf), sizeof(pointsToLeaf));

    // leer cantidad de hijos
    size_t numChildren;
    in.read(reinterpret_cast<char *>(&numChildren), sizeof(numChildren));

    // leer hijos
    for (size_t i = 0; i < numChildren; ++i) {
        SsNode *child = pointsToLeaf ? static_cast<SsNode *>(new SsLeaf(D)) : static_cast<SsNode *>(new SsInnerNode(D));
        child->loadFromStream(in, this);
        children.push_back(child);
    }
}

void SsLeaf::loadFromStream(std::istream &in, SsNode *parent) {
    this->parent = parent;

    // Leer centroid
    centroid.readFromFile(in, D);

    // Leer radio
    float radius_ = 0;
    in.read(reinterpret_cast<char *>(&radius_), sizeof(radius_));
    this->radius = radius_;

    // Leer numero de puntos
    size_t numPoints;
    in.read(reinterpret_cast<char *>(&numPoints), sizeof(numPoints));

    // Leer puntos
    points.resize(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
        points[i].readFromFile(in, D);
    }

    // Leer rutas (paths)
    size_t numPaths;
    in.read(reinterpret_cast<char *>(&numPaths), sizeof(numPaths));
    paths.resize(numPaths);
    for (size_t i = 0; i < numPaths; ++i) {
        size_t pathLength;
        in.read(reinterpret_cast<char *>(&pathLength), sizeof(pathLength));
        char *buffer = new char[pathLength + 1];
        in.read(buffer, (long) pathLength);
        buffer[pathLength] = '\0';
        paths[i] = std::string(buffer);
        delete[] buffer;
    }
}

void SsTree::saveToFile(const std::string &filename) const {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Cannot open file for writing");
    }

    // Guardar las dimensiones de la estructura
    out.write(reinterpret_cast<const char *>(&D), sizeof(D));

    // Guardar si el root es hija o nodo interno
    bool isLeaf = root->isLeaf();
    out.write(reinterpret_cast<const char *>(&isLeaf), sizeof(isLeaf));

    // Guardar el resto de la estructura
    root->saveToStream(out);
    out.close();
}

void SsTree::loadFromFile(const std::string &filename) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot open file for reading");
    }
    if (root) {
        delete root;
        root = nullptr;
    }

    // Aquí se asume que el primer valor determina las dimensiones
    in.read(reinterpret_cast<char *>(&D), sizeof(D));

    // El segundo valor determina si el root es hoja
    bool isLeaf;
    in.read(reinterpret_cast<char *>(&isLeaf), sizeof(isLeaf));
    if (isLeaf) {
        root = new SsLeaf(D);
    } else {
        root = new SsInnerNode(D);
    }
    root->loadFromStream(in, nullptr);
    in.close();
}

std::vector<Point> SsInnerNode::getEntriesCentroids() const {
    std::vector<Point> centroids;
    centroids.reserve(children.size());

    for (SsNode *node: children) {
        centroids.push_back(node->centroid);
    }

    return centroids;
}

void SsInnerNode::sortEntriesByCoordinate(size_t coordinateIndex) {
    std::sort(children.begin(), children.end(), [&](SsNode *n1, SsNode *n2) {
        return n1->centroid[coordinateIndex] < n2->centroid[coordinateIndex];
    });
}

SsInnerNode::SsInnerNode(std::size_t dims) : SsNode(dims) {
}

std::pair<SsNode *, SsNode *> SsInnerNode::split() {
    long splitIndex = findSplitIndex();

    SsNode *newNode1 = new SsInnerNode(D);
    dynamic_cast<SsInnerNode *>(newNode1)->children = std::vector<SsNode *>(this->children.begin(),
                                                                            children.begin() + splitIndex);

    SsNode *newNode2 = new SsInnerNode(D);
    dynamic_cast<SsInnerNode *>(newNode2)->children = std::vector<SsNode *>(this->children.begin() + splitIndex,
                                                                            this->children.end());

    newNode1->updateBoundingEnvelope();
    newNode2->updateBoundingEnvelope();

    newNode1->parent = this->parent;
    newNode2->parent = this->parent;

    return std::make_pair(newNode1, newNode2);
}

SsNode *SsInnerNode::findMinRadiusIncreaseChild(const Point &target) const {
    SsNode *closest = nullptr;
    NType min_radius_increase = NType::max_value();
    NType min_radius = NType::max_value();

    for (SsNode *node: children) {
        NType dis = distance(node->centroid, target);

        if (dis <= node->radius) {
            if (min_radius > (float) node->getEntriesSize()) {
                min_radius = (float) node->getEntriesSize();
                closest = node;
            }

            min_radius_increase = 0;
        } else {
            if (min_radius_increase == 0) {
                continue;
            }

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

bool SsInnerNode::isLeaf() const {
    return false;
}

std::vector<Point> SsInnerNode::buildReduction() {
    std::vector<Point> result;
    result.reserve(2 * D * children.size());

    for (SsNode *node: this->children) {
        for (int i = 0; i < D; ++i) {
            Point p = node->centroid;
            p[i] += node->radius;
            result.push_back(p);
            p[i] -= node->radius * 2.0;
            result.push_back(p);
        }
    }

    return result;
}

void SsInnerNode::updateBoundingEnvelope() {
    std::vector<Point> centroids = getEntriesCentroids();
    std::vector<Point> reduction = buildReduction();

    Point mean = meanPoint(centroids);

    auto frankWolfeAlgorithm = FrankWolfeAlgorithm(reduction, mean, frankWolfeIterations);
    for (SsNode *node: children) {
        NType dis = distance(frankWolfeAlgorithm.second, node->centroid) + node->radius;
        if (dis > frankWolfeAlgorithm.first) {
            frankWolfeAlgorithm.first = dis;
        }
    }

    this->radius = frankWolfeAlgorithm.first;
    this->centroid = frankWolfeAlgorithm.second;
}

std::pair<SsNode *, SsNode *> SsInnerNode::insert(const Point &point, const std::string &path) {
    SsNode *closestChild = findMinRadiusIncreaseChild(point);

    NType prev_radius = closestChild->radius;
    std::pair<SsNode *, SsNode *> newChildren = closestChild->insert(point, path);
    SsNode *newChild1 = newChildren.first;
    SsNode *newChild2 = newChildren.second;

    if (newChild1 == nullptr) {
        if (closestChild->radius > prev_radius) {
            children.erase(std::remove(children.begin(), children.end(), closestChild),
                           children.end());
            children.emplace_back(closestChild);
            updateBoundingEnvelope();
        }
        return std::make_pair(nullptr, nullptr);
    } else {
        children.erase(std::remove(children.begin(), children.end(), closestChild),
                       children.end());
        children.emplace_back(newChild1);
        children.emplace_back(newChild2);
        updateBoundingEnvelope();

        if (children.size() <= Settings::M) {
            return std::make_pair(nullptr, nullptr);
        }
    }

    return split();
}

int SsInnerNode::getEntriesSize() const {
    return (int) this->children.size();
}

int SsLeaf::getEntriesSize() const {
    return (int) this->points.size();
}

std::vector<Point> SsLeaf::getEntriesCentroids() const {
    return points;
}

void SsLeaf::sortEntriesByCoordinate(size_t coordinateIndex) {
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

SsLeaf::SsLeaf(size_t dims) : SsNode(dims) {}

std::pair<SsNode *, SsNode *> SsLeaf::split() {
    size_t splitIndex = findSplitIndex();

    auto *newNode1 = new SsLeaf(D);
    newNode1->points.reserve(splitIndex);
    for (size_t i = 0; i < splitIndex; ++i) {
        newNode1->points.push_back(points[i]);
        newNode1->paths.push_back(paths[i]);
    }
    newNode1->parent = parent;

    auto *newNode2 = new SsLeaf(D);
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

bool SsLeaf::isLeaf() const {
    return true;
}

void SsLeaf::updateBoundingEnvelope() {
    if (points.size() == 1) {
        this->radius = 0;
        this->centroid = points[0];
        return;
    } else if (points.size() == 2) {
        this->centroid = (points[0] + points[1]) / 2.0;
        this->radius = distance(points[0], centroid);
        return;
    }

    Point mean = meanPoint(points);
    auto frankWolfeAlgorithm = FrankWolfeAlgorithm(points, mean, frankWolfeIterations);

    this->radius = frankWolfeAlgorithm.first;
    this->centroid = frankWolfeAlgorithm.second;
}

std::pair<SsNode *, SsNode *> SsLeaf::insert(const Point &point, const std::string &path) {
    if (std::find(points.begin(), points.end(), point) != points.end()) {
        return std::make_pair(nullptr, nullptr);
    }

    points.push_back(point);
    paths.push_back(path);

    if (distance(point, centroid) > radius) {
        updateBoundingEnvelope();
    }

    if (points.size() <= Settings::M) {
        return std::make_pair(nullptr, nullptr);
    }

    return split();
}


void SsTree::kNNQuery(const Point &target, size_t k, SsNode *node, NType radius,
                       std::priority_queue<std::pair<NType, std::pair<Point, std::string>>, std::vector<std::pair<NType, std::pair<Point, std::string>>>, decltype(compare)> &result) const {
    if (node->isLeaf()) {
        auto *leafNode = dynamic_cast<SsLeaf *>(node);
        for (int i = 0; i < leafNode->points.size(); ++i) {
            Point p = leafNode->points[i];

            NType dist = std::max(NType(0), distance(target, p));
            // Regla 2
            if (radius + std::max(NType(0), dist) < std::max(NType(0), distance(target, leafNode->centroid))) {
                continue;
            }
            // Regla 4
            if (radius + std::max(NType(0), distance(target, leafNode->centroid)) <
                std::max(NType(0), distance(p, leafNode->centroid))) {
                continue;
            }

            if (dist < radius) {
                result.push({dist, {p, leafNode->paths[i]}});
                if (result.size() > k) {
                    result.pop();
                }
            }
        }
    } else {
        auto *innerNode = dynamic_cast<SsInnerNode *>(node);
        // Se obtiene la menor distancia hacia un nodo
        for (auto *child: innerNode->children) {
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
                kNNQuery(target, k, child, radius, result);
            }
        }
    }
}

SsTree::SsTree(std::size_t dims) : root(nullptr), D(dims) {
}

SsTree::~SsTree() {
    delete root;
}

void SsTree::insert(const Point &point, const std::string &path) {
    if (root == nullptr) {
        root = new SsLeaf(D);
        root->parent = nullptr;
        root->centroid = point;
        root->radius = 0;
    }

    std::pair<SsNode *, SsNode *> newChildren = root->insert(point, path);
    SsNode *newChild1 = newChildren.first;
    SsNode *newChild2 = newChildren.second;

    if (newChild1 == nullptr) {
        return;
    }

    root = new SsInnerNode(D);
    root->parent = nullptr;
    newChild1->parent = root;
    newChild2->parent = root;

    dynamic_cast<SsInnerNode *>(root)->children.emplace_back(newChild1);
    dynamic_cast<SsInnerNode *>(root)->children.emplace_back(newChild2);
    (dynamic_cast<SsInnerNode *>(root))->updateBoundingEnvelope();
}

std::vector<std::string> SsTree::kNNQuery(const Point &center, size_t k) const {
    if (root == nullptr) return {};
    std::priority_queue<std::pair<NType, std::pair<Point, std::string>>, std::vector<std::pair<NType, std::pair<Point, std::string>>>, decltype(compare)> result{
            compare};
    this->kNNQuery(center, k, this->root, NType::max_value(), result);
    std::vector<std::string> points;
    for (int i = 0; i < k; ++i) {
        points.push_back(result.top().second.second);
        result.pop();
    }
    std::vector<std::string> points_inv;
    std::reverse_copy(points.begin(), points.end(), std::back_inserter(points_inv));
    return points_inv;
}

std::pair<NType, Point> FrankWolfeAlgorithm(std::vector<Point> &points, Point &seed, int k) {
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

Point meanPoint(const std::vector<Point> &points) {
    Point mean = points[0];

    for (int i = 1; i < points.size(); ++i) {
        mean += points[i];
    }

    mean /= (float) points.size();
    return mean;
}
