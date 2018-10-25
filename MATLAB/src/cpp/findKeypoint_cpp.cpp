

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>
//#include <eigen>



using namespace matlab::data;
using matlab::mex::ArgumentList;



class MexFunction : public matlab::mex::Function {
public:

  void operator()(ArgumentList outputs, ArgumentList inputs) {
    checkArguments(outputs, inputs);

    matlab::data::TypedArray<double> unit = std::move(inputs[0]);
    outputs[0] = unit;

    matlab::data::CellArray NeighborsOut = std::move(inputs[1]);
    outputs[1] = NeighborsOut;


    // int a = GetDimensions(NeighborsOut);
    // matlab::data::TypedArrayRef<double> ref NeighborsOut[0];

  }







  void checkArguments(ArgumentList outputs, ArgumentList inputs){

    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    matlab::data::ArrayFactory factory;

    // Check first input argument is a matrix (DoG)
    if (inputs[0].getType() != matlab::data::ArrayType::DOUBLE ||
        inputs[0].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) {
          matlabPtr->feval(u"error", 0,
                           std::vector<matlab::data::Array>({ factory.createScalar("Input must be double array") }));
        }

    // Check first input argument is a matrix (Neighbors)
    if (inputs[1].getType() != matlab::data::ArrayType::CELL) {
          matlabPtr->feval(u"error", 0,
                           std::vector<matlab::data::Array>({ factory.createScalar("Input must be cell array") }));
        }



    // Check number of output ArgumentList

    if (outputs.size() > 2) {
      matlabPtr->feval(u"error", 0,
                       std::vector<matlab::data::Array>( {factory.createScalar("Only two outputs are returned.")} ));
    }
  }
};
