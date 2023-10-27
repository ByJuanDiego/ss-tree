#include <iostream>
#include <vector>
#include "SStree.h"
#include <H5Cpp.h>

struct ImageData {
    Point embedding;
    std::string path;

    friend std::ostream& operator << (std::ostream& os, const ImageData& img) {
        os << img.path << ": " << img.embedding;
        return os;
    }
};

std::vector<ImageData> readEmbeddingsFromHDF5(const H5std_string &FILE_NAME) {
    std::vector<ImageData> data;
    try {
        H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);
        H5::DataSet featuresDataset = file.openDataSet("features");
        H5::DataSpace featuresSpace = featuresDataset.getSpace();
        const int num_embeddings = featuresSpace.getSimpleExtentNpoints();
        auto *embeddingData = new double[num_embeddings];
        featuresDataset.read(embeddingData, H5::PredType::NATIVE_DOUBLE);
        auto itemType = H5::PredType::NATIVE_DOUBLE;
        auto memType = H5::VarLenType(&itemType);
        H5::DataSet dataset = file.openDataSet("paths");
        H5::DataSpace dataSpace = dataset.getSpace();
        hsize_t rank;
        hsize_t dims[1];
        rank = dataSpace.getSimpleExtentDims(dims);
        hsize_t memDims[1] = {1};
        H5::DataSpace memspace(rank, memDims);
        hsize_t dataCount[1];
        hsize_t dataOffset[1];
        hsize_t memCount[1];
        hsize_t memOffset[1];
        hsize_t pointDims = num_embeddings / dims[0];
        for (hsize_t i = 0; i < dims[0]; i++) {
            dataCount[0] = 1;
            dataOffset[0] = i;
            memCount[0] = 1;
            memOffset[0] = 0;
            dataSpace.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);
            memspace.selectHyperslab(H5S_SELECT_SET, memCount, memOffset);
            auto *rdata = new hvl_t[1];
            dataset.read(rdata, memType, memspace, dataSpace);
            auto *ptr = (double *) rdata[0].p;
            std::string path;
            // Leer path
            for (int j = 0; j < rdata[0].len; j++) {
                auto *val = (double *) &ptr[j];
                path.push_back((char) *val);
            }
            // Leer embedding
            Point point(pointDims);
            for (int j = 0; j < pointDims; ++j) {
                point[j] = (float) embeddingData[i * pointDims + j];
            }
            data.push_back({point, path});
            delete [] rdata;
        }
        delete [] embeddingData;
        file.close();
    } catch (H5::Exception &error) {
        std::cerr << error.getCDetailMsg() << std::endl;
    }
    return data;
}

int main() {
    const H5std_string FILE_NAME("embedding.hdf5");
    std::vector<ImageData> data = readEmbeddingsFromHDF5(FILE_NAME);

    SsTree tree(448);
//    tree.loadFromFile("tree4.dat");
//    std::cout << tree.count() << std::endl;
//    tree.test();
//    tree.print();

    int i = 0;
    for (const ImageData& img: data) {
        tree.insert(img.embedding, img.path);
        ++i;
        if (!(i % 1000)) {
            std::cout << i << ": ";
            tree.print();
        }
    }

    tree.print();
    tree.saveToFile("tree5.dat");


    return 0;
}

// tree.dat: Contiene el metodo de mean-shift y kfn para leaf nodes, para inner nodes tiene el metodo de MBB, MBB2, IMBB y mean point sobre la direccion de su media
// tree2.dat: Contiene el metodo de mean-shift y kfn para leaf nodes & el de desplazamiento desde el mean point y el mean point sin considerar el nuevo punto hacia el nuevo punto, para inner nodes tiene el metodo de MBB, MBB2, IMBB y mean point sobre la direccion de su media
// los principales cambios respecto al 2 son los cambios propuestos por Remon... "avanzar" en direccion de la cosa insertada
// nueva logica en update bounding envelope para splits, calcular el punto medio entre los dos circulos tomando en cuenta su radio
// desempate por el menor candidad de entries en el findminradius (no sirvio), se descarta
// desempatar por el de menor radio en caso todos engloben al punto

