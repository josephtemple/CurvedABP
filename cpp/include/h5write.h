// to output to HDF5 file format
#pragma once
#include <H5Cpp.h>

// Template method for writing to HDF5
template<typename T>
struct H5TypeMap;

template<>
struct H5TypeMap<double> {
    static H5::PredType type() { return H5::PredType::NATIVE_DOUBLE; }
};

template<>
struct H5TypeMap<unsigned int> {
    static H5::PredType type() { return H5::PredType::NATIVE_UINT; }
};

template<>
struct H5TypeMap<int> {
    static H5::PredType type() { return H5::PredType::NATIVE_INT; }
};

// template method for writing to a given container
template<typename T, typename H5Container>
void write_H5(const std::string& name, H5Container& write_to, const T& val) {
    hsize_t dims[1] = {1};
    H5::DataSpace ds(1, dims);

    auto dtype = H5TypeMap<T>::type();

    H5::DataSet dset = write_to.createDataSet(name, dtype, ds);
    dset.write(&val, dtype);
}