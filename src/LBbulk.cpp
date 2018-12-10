#include "LBbulk.h"


Bulk::Bulk(int nBulkNodes)
{
  nBulkNodes_=nBulkNodes;
  bulkNode_ = new int[nBulkNodes_];
}

Bulk::~Bulk()
{
  delete [] bulkNode_;
}
