#ifndef LBBULK_H
#define LBBULK_H


#include <iostream>
#include "LBglobal.h"
#include "LBfield.h"
#include "LBgrid.h"
#include "LBd2q9.h"
#include "LBboundary.h"

class Bulk
{
public:
  Bulk(int nBulkNodes);
  ~Bulk();

protected:
  int nBulkNodes_;
  int* bulkNode_;
  
};






#endif // LBBOUNDARY_H
