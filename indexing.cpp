#include <iostream>
#include <vector>
#include "SStree.h"
#include <nlohmann/json.hpp>
#include <fstream>

struct ImageData {
    std::vector<Point> embeddings;
    std::vector<std::string> paths;

};

ImageData readEmbeddingsFromJson(const std::string& FILE_NAME) {
    ImageData data;
    data.embeddings.reserve(26179);
    data.paths.reserve(26179);

    try {
        std::ifstream file(FILE_NAME);
        if (!file.is_open()) {
            throw std::runtime_error("Unable to open JSON file.");
        }

        nlohmann::json jsonData;
        file >> jsonData;
        for (auto it = jsonData.begin(); it != jsonData.end(); ++it) {
            const std::string& path = it.key();
            std::vector<float> values = it.value();
            std::vector<NType> embedding;
            embedding.reserve(values.size());
            for (auto& x: values) {
                embedding.emplace_back(x);
            }

            data.embeddings.emplace_back(embedding);
            data.paths.push_back(path);
        }
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return data;
}

int main() {
    const std::string FILE_NAME("embedding.json");
    ImageData data = readEmbeddingsFromJson(FILE_NAME);
    SsTree tree(448);

    std::cout << "Begin\n";

    for (int i = 0; i < data.embeddings.size(); ++i) {
        tree.insert(data.embeddings[i], data.paths[i]);

        if ((i+1)%100 == 0) {
            std::cout << (i+1) << ": " << "\n";
            tree.print();
            tree.test();
        }
    }

    std::cout << "End\n";
    tree.print();
    tree.test();

    tree.saveToFile("embedding1.dat");
}
