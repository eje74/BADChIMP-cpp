#include "LBbulk.h"


Bulk::Bulk(const int nBulkNodes) : nBulkNodes_(nBulkNodes), bulkNode_(nBulkNodes_)
{
  //nBulkNodes_=nBulkNodes;
  nAddedNodes_ = 0;
  //bulkNode_ = new int[nBulkNodes_];
}

/*
Bulk::~Bulk()
{
  delete [] bulkNode_;
}
*/

void Bulk::addBulkNode(const int nodeNo )
{
  if (nAddedNodes_ < nBulkNodes_) {
        bulkNode_[nAddedNodes_] = nodeNo;
    }
    else {
        std::cout << "ERROR: Added to many boundary nodes!" << std::endl;
        exit(1);
    }
  nAddedNodes_ += 1;
}

