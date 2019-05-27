#ifndef LBBULK_H
#define LBBULK_H


#include <iostream>
#include "LBglobal.h"
#include "LBfield.h"
#include "LBgrid.h"
//#include "LBd2q9.h"
#include "LBlatticetypes.h"
#include "LBboundary.h"

class Bulk
{
public:
  Bulk(int nBulkNodes);
  ~Bulk();

  void addBulkNode( int nodeNo );
  int nElements() const;
  int nodeNo(const int elementNo) const;

private:
  int nBulkNodes_;
  int nAddedNodes_;// Counter for number of added boundary nodes
  //int* bulkNode_;
  std::vector<int> bulkNode_;
  
};

// Getter for number of elements
inline int Bulk::nElements() const
{
    return nBulkNodes_;
}

// Returns the node number from bulk element (without ghost nodes)
inline int Bulk::nodeNo(const int elementNo) const
{
    return bulkNode_[elementNo];
}





#endif // LBBOUNDARY_H
