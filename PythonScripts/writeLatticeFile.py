# -*- coding: utf-8 -*-
# WRITE AN LATTICE CLASS USING OLD BADCHIMP LATTICE INPUT FILES


def int_x_vec(i, v_name, v_index):
    if i == 0:
        return ""
    if i == 1:
        return " +" + v_name + "[" + str(v_index) + "]"
    if i == -1:
        return " -" + v_name + "[" + str(v_index) + "]"
    return str(i) + "*" + v_name + "[" + str(v_index) + "]"


def write_dot(dxqy, nd):
    cl = "inline lbBase_t {0:s}::dot(const lbBase_t* leftVec, const lbBase_t* rightVec)".format(dxqy)
    print(cl)
    f.write(cl+"\n")
    cl = "{"
    print(cl)
    f.write(cl+"\n")
    cl = "    return "
    for d in range(nd):
        cl += "leftVec[{0:d}]*rightVec[{0:d}]".format(d)
        if d != nd-1:
            cl += " + "
        else:
            cl += ";"
    print(cl)
    f.write(cl+"\n")
    cl = "}"
    print(cl)
    f.write(cl+"\n"+"\n")
    print("")

# inline lbBase_t D2Q9::cDot(const int qDir, const lbBase_t* rightVec)
# {
#     return c(qDir, 0)*rightVec[0] + c(qDir, 1)*rightVec[1];
# }
def write_cDot(dxqy, nd):
    cl = "inline lbBase_t {0:s}::cDot(const int qDir, const lbBase_t* rightVec)".format(dxqy)
    print(cl)
    f.write(cl+"\n")
    cl = "{"
    print(cl)
    f.write(cl+"\n")
    cl = "    return "
    for d in range(nd):
        cl += "c(qDir, {0:d})*rightVec[{0:d}]".format(d)
        if d != nd-1:
            cl += " + "
        else:
            cl += ";"
    print(cl)
    f.write(cl+"\n")
    cl = "}"
    print(cl)
    f.write(cl+"\n"+"\n")
    print("")

# inline void D2Q9::cDotAll(const lbBase_t* vec, lbBase_t* ret)
# {
#     ret[0] = vec[0];
#     ret[1] = vec[0] + vec[1];
#     ret[2] = vec[1];
#     ret[3] = vec[1] - vec[0];
#     ret[4] = -vec[0];
#     ret[5] = -vec[0] - vec[1];
#     ret[6] = -vec[1];
#     ret[7] = vec[0] - vec[1];
#     ret[8] = 0;
# }
#
def write_cDotAll(dxqy, nd, nq, cv):
    cl = "inline void {0:s}::cDotAll(const lbBase_t* vec, lbBase_t* ret)".format(dxqy)
    print(cl)
    f.write(cl+"\n")
    cl = "{"
    print(cl)
    f.write(cl+"\n")
    for q in range(nq):
        cl = "ret[{0:d}] =".format(q)
        if all([x == 0 for x in cv[q]]):
            cl += " 0.0;"
        else:
            for d in range(nd):
                cl += int_x_vec(cv[q][d], "vec", d)
            cl += ";"
        print(cl)
        f.write(cl+"\n")
    cl = "}"
    print(cl)
    f.write(cl+"\n"+"\n")
    print("")
    
# inline void D2Q9::grad(const lbBase_t* rho, lbBase_t* ret)
# {
#     ret[0] =  w1c2Inv * (rho[0] - rho[4]);
#     ret[0] += w2c2Inv * (rho[1] - rho[3] - rho[5] + rho[7]);
#     ret[1] =  w1c2Inv * (rho[2] - rho[6]);
#     ret[1] += w2c2Inv * (rho[1] + rho[3] - rho[5] - rho[7]);
# }
#
def write_grad(dxqy, nd, nq, cv, cL):
    cl = "inline void {0:s}::grad(const lbBase_t* rho, lbBase_t* ret)".format(dxqy)
    print(cl)
    f.write(cl+"\n")
    cl = "{"
    print(cl)
    f.write(cl+"\n")
    for d in range(nd):
        cl = "ret[{0:d}] =".format(d)
        for clength in range(max(cL) + 1):
            qval = [q for q in range(nq) if cv[q][d] != 0 and cL[q] == clength]
            if len(qval) > 0:
                cl += "+ w{0:d}c2Inv * (".format(clength)
                for qv in qval:
                    cl += int_x_vec(cv[qv][d], " rho", qv)
                cl += " ) "
        cl += ";"
        print(cl)
        f.write(cl+"\n")
    cl = "}"
    print(cl)
    f.write(cl+"\n"+"\n")
    print("")
    
# inline void D2Q9::qSum(const lbBase_t* dist, lbBase_t& ret)
# {
#     ret = 0.0;
#     for (int q = 0; q < nQ; ++q)
#         ret += dist[q];
# }

def write_qSum(dxqy):
    cl = "inline void {0:s}::qSum(const lbBase_t* dist, lbBase_t& ret)".format(dxqy)
    print(cl)
    f.write(cl+"\n")
    cl = "{"
    print(cl)
    f.write(cl+"\n")
    cl = "ret = 0.0;"    
    print(cl)
    f.write(cl+"\n")
    cl = "for (int q = 0; q < nQ; ++q)"
    print(cl)
    f.write(cl+"\n")
    cl = "ret += dist[q];"
    print(cl)
    f.write(cl+"\n")
    cl = "}"
    print(cl)
    f.write(cl+"\n"+"\n")
    print("")
    
# inline void D2Q9::qSumC(const lbBase_t* dist, lbBase_t* ret)
# {
#     ret[0] = dist[0] + dist[1]            - dist[3] - dist[4]  - dist[5]           + dist[7];
#     ret[1] =           dist[1] + dist[2]  + dist[3]            - dist[5] - dist[6] - dist[7];
# }

def write_qSumC(dxqy, nd, nq, cv):
    cl = "inline void {0:s}::qSumC(const lbBase_t* dist, lbBase_t* ret)".format(dxqy)
    print(cl)
    f.write(cl+"\n")
    cl = "{"
    print(cl)
    f.write(cl+"\n")
    for d in range(nd):
        cl = "ret[{0:d}] =".format(d)
        for q in range(nq):
            if any([x != 0 for x in cv[q]]):
                if(cv[q][d]!=0):
                    cl += int_x_vec(cv[q][d]," dist",q)
        print(cl + ";")
        f.write(cl + ";" + "\n")
    cl = "}"
    print(cl)
    f.write(cl+"\n"+"\n")
    print("")
    
    
# inline void D2Q9::gradPush(const lbBase_t &scalarVal, const int *neigList, VectorField<D2Q9> &grad)
# {
#     const lbBase_t valTmp1  = scalarVal * c2Inv * w1;
#     const lbBase_t valTmp2  = scalarVal * c2Inv * w2;
#
#     int nodeNeigNo = neigList[0];
#     grad(0,0,nodeNeigNo) -= valTmp1;
#
#     nodeNeigNo = neigList[1];
#     grad(0,0,nodeNeigNo) -= valTmp2;
#     grad(0,1,nodeNeigNo) -= valTmp2;
#
#     nodeNeigNo = neigList[2];
#     grad(0,1,nodeNeigNo) -= valTmp1;
#
#     nodeNeigNo = neigList[3];
#     grad(0,0,nodeNeigNo) += valTmp2;
#     grad(0,1,nodeNeigNo) -= valTmp2;
#
#     nodeNeigNo = neigList[4];
#     grad(0,0,nodeNeigNo) += valTmp1;
#
#     nodeNeigNo = neigList[5];
#     grad(0,0,nodeNeigNo) += valTmp2;
#     grad(0,1,nodeNeigNo) += valTmp2;
#
#     nodeNeigNo = neigList[6];
#     grad(0,1,nodeNeigNo) += valTmp1;
#
#     nodeNeigNo = neigList[7];
#     grad(0,0,nodeNeigNo) -= valTmp2;
#     grad(0,1,nodeNeigNo) += valTmp2;
# }

def write_gradPush(dxqy, nd, nq, cv, cL):
    cl = "inline void {0:s}::gradPush(const lbBase_t &scalarVal, const int *neigList, VectorField<{0:s}> &grad)".format(dxqy)
    print(cl)
    f.write(cl+"\n")
    cl = "{"
    print(cl)
    f.write(cl+"\n")
    
    for clength in range(max(cL) + 1):
        if clength > 0:
            cl = "const lbBase_t valTmp{0:d}  = scalarVal * c2Inv * w{0:d};".format(clength)
        #cl = "const lbBase_t valTmp{0:d}  = scalarVal * c2Inv * w{0:d};".format(x)
            print(cl)
            f.write(cl+"\n")  
    print("")   
    f.write("\n")
    
    for q in range(nq-1):
        cl = ""
        if q == 0:
            cl += "int "
        cl += "nodeNeigNo = neigList[{0:d}];".format(q)
        print(cl)
        f.write(cl+"\n")
        for d in range(nd):
            if cv[q][d]!=0:
                cl = "grad(0,{0:d},nodeNeigNo) ".format(d)
                if cv[q][d]>0:
                    cl += "-= "
                else:
                    cl += "+= "
                cl += "valTmp{0:d};".format(cL[q])
                print(cl)
                f.write(cl+"\n")
        print("") 
        f.write("\n")   
    cl = "}"
    print(cl)
    f.write(cl+"\n"+"\n")
    print("")
    
    
#------------------------------------------------------------------------
#----------------------------LBdXqY.cpp----------------------------------
#------------------------------------------------------------------------

#Number of dimensions
nD = 2
#Number of lattice directions
nQ = 9

#Number of lattice pairs
nDirPairs = 4
#Number of lattice directions pointing to neighbors
nQNonZero = 8

# 1st lattice constant
c2Inv = 3.0
# 2nd lattice constant
c4Inv = 9.0

# weight fractions:
# weight denominator
wGcd = 36
#weight numerator
#Order: rest, lenght 1, lenght 2,...
wN = [16, 4, 1]

#Lattice vectors:
vec = []
#x-component
vecDimX = [1, 1, 0, -1, -1, -1, 0, 1, 0]
vec.append(vecDimX)
#y-component
vecDimY = [0, 1, 1, 1, 0, -1, -1, -1, 0]
vec.append(vecDimY)


# // Two phase values
# B weights used in surface tension
#fractions:
# denominator
bGcd = 108
#numerator
bN = [-16, 8, 5]

#-------------------------------------------------------------

#Number of dimensions
nD = 3
#Number of lattice directions
nQ = 19

#Number of lattice pairs
nDirPairs = 9
#Number of lattice directions pointing to neighbors
nQNonZero = 18

# 1st lattice constant
c2Inv = 3.0
# 2nd lattice constant
c4Inv = 9.0

# weight fractions:
# weight denominator
wGcd = 36
#weight numerator
#Order: rest, lenght 1, lenght 2,...
wN = [12, 2, 1]

#Lattice vectors:
vec = []
#x-component
vecDimX = [1, 0, 0, 1, 1, 1, 1, 0, 0, -1, 0, 0, -1, -1, -1, -1, 0, 0, 0]
vec.append(vecDimX)
#y-component
vecDimY = [0, 1, 0, 1, -1, 0, 0, 1, 1, 0, -1, 0, -1, 1, 0, 0, -1, -1, 0]
vec.append(vecDimY)
#z-component
vecDimZ = [0, 0, 1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 0, 0, -1, 1, -1, 1, 0]
vec.append(vecDimZ)


# // Two phase values
# B weights used in surface tension
#fractions:
# denominator
bGcd = 54
#numerator
bN = [-12, 1, 2]

f=open("LBd{0:d}q{1:d}.h".format(nD, nQ),"w+")

latticeName = "D{0:d}Q{1:d}".format(nD,nQ)

codeLine = "#ifndef LBD{0:d}Q{1:d}_H".format(nD,nQ)
print(codeLine)
f.write(codeLine+"\n")
codeLine = "#define LBD{0:d}Q{1:d}_H".format(nD,nQ)
print(codeLine)
f.write(codeLine+"\n")
f.write("\n")

codeLine = '#include "LBglobal.h"'
print(codeLine)
f.write(codeLine+"\n")
codeLine = '#include "LBfield.h"'
print(codeLine)
f.write(codeLine+"\n")
f.write("\n")
f.write("\n")

codeLine = '// See "LBlatticetypes.h" for description of the structure'
print(codeLine)
f.write(codeLine+"\n")
f.write("\n")

codeLine = "struct " + latticeName + "{"
print(codeLine)
f.write(codeLine+"\n")
f.write("\n")

codeLine = "static constexpr int nD = " + str(nD) + ";"
print(codeLine)
f.write(codeLine+"\n")
codeLine = "static constexpr int nQ = " + str(nQ) + ";"
print(codeLine)
f.write(codeLine+"\n")

codeLine = "static constexpr int nDirPairs_ = " + str(nDirPairs) + ";"
print(codeLine)
f.write(codeLine+"\n")
codeLine = "static constexpr int nQNonZero_ = " + str(nQNonZero) + ";"
print(codeLine)
f.write(codeLine+"\n")
f.write("\n")


codeLine = "static constexpr lbBase_t c2Inv = " + str(c2Inv) + ";"
print(codeLine)
f.write(codeLine+"\n")
codeLine = "static constexpr lbBase_t c4Inv = " + str(c4Inv) + ";"
print(codeLine)
f.write(codeLine+"\n")
codeLine = "static constexpr lbBase_t c2 = 1.0 / c2Inv;"
print(codeLine)
f.write(codeLine+"\n")
codeLine = "static constexpr lbBase_t c4 = 1.0 / c4Inv;"
print(codeLine)
f.write(codeLine+"\n")
codeLine = "static constexpr lbBase_t c4Inv0_5 = 0.5 * c4Inv;"
print(codeLine)
f.write(codeLine+"\n")
f.write("\n")


for num, wn in enumerate(wN):
    codeLine = "static constexpr lbBase_t w{0:d} = {1:.1f}/{2:.1f};".format(num, wn, wGcd)
    print(codeLine)
    f.write(codeLine+"\n")
    codeLine = "static constexpr lbBase_t w{0:d}c2Inv = w{0:d}*c2Inv;".format(num)
    print(codeLine)
    f.write(codeLine+"\n")
f.write("\n")

#### vector
# // Remember do define the static arrays in the cpp file as well.


cBasis = []
cLength = []
for q in range(nQ):
    cvec = []
    cSquare = 0
    for d in range(nD):
        cSquare += vec[d][q]**2
        cvec.append(vec[d][q])
    cLength.append(cSquare)
    cBasis.append(cvec)


codeLineW = "static constexpr lbBase_t w[{0:d}] = ".format(nQ) + "{"
codeLineCVec = "static constexpr int cDMajor_[{0:d}] = ".format(nQ*nD) + "{"
codeLineCNorm = "static constexpr lbBase_t cNorm[{0:d}] = ".format(nQ) + "{"
for q in range(nQ):
    codeLineW += "w{0:d}".format(cLength[q])
    if cLength[q] == 0:
        codeLineCNorm += "0.0"
    elif cLength[q] == 1:
        codeLineCNorm += "1.0"
    else:
        codeLineCNorm += "SQRT{0:d}".format(cLength[q])
    for d in range(nD-1):
        codeLineCVec += "{0:d}, ".format(vec[d][q])
    codeLineCVec += "{0:d}".format(vec[nD-1][q])
    if q != nQ-1:
        codeLineW += ", "
        codeLineCVec += ", "
        codeLineCNorm += ", "
    else:
        codeLineW += "};"
        codeLineCVec += "};"
        codeLineCNorm += "};"
print(codeLineW)
f.write(codeLineW+"\n")
print(codeLineCVec)
f.write(codeLineCVec+"\n")
print(codeLineCNorm)
f.write(codeLineCNorm+"\n")

# // Two phase values




codeLine = ""
for num, bn in enumerate(bN):
    codeLine = "static constexpr lbBase_t B{0:d} = {1:.1f}/{2:.1f};".format(num, bn, bGcd)
    print(codeLine)
    f.write(codeLine+"\n")

codeLine = "static constexpr lbBase_t B[{0:d}] = ".format(nQ) + "{"
for q in range(nQ):
    codeLine += "B{0:d}".format(cLength[q])
    if q != nQ-1:
        codeLine += ", "
    else:
        codeLine += "};"

print(codeLine)
f.write(codeLine+"\n"+"\n")
print("")
codeLine = "// Functions"
print(codeLine)
f.write(codeLine+"\n"+"\n")
print("")
codeLine = "inline static int c(const int qDirection, const int dimension)  {return cDMajor_[nD*qDirection + dimension];}"
print(codeLine)
f.write(codeLine+"\n")
codeLine = "inline static int reverseDirection(const int qDirection) {return (qDirection + nDirPairs_) % nQNonZero_;}"
print(codeLine)
f.write(codeLine+"\n")
f.write("\n")

codeLine = "static lbBase_t dot(const lbBase_t* leftVec, const lbBase_t* rightVec);"
print(codeLine)
f.write(codeLine+"\n")
codeLine = "static lbBase_t cDot(const int qDir, const lbBase_t* rightVec);"
print(codeLine)
f.write(codeLine+"\n")
codeLine = "static void cDotAll(const lbBase_t* vec, lbBase_t* ret);"
print(codeLine)
f.write(codeLine+"\n")
codeLine = "static void grad(const lbBase_t* rho, lbBase_t* ret);"
print(codeLine)
f.write(codeLine+"\n")
f.write("\n")
codeLine = "static void qSum(const lbBase_t* dist, lbBase_t& ret);"
print(codeLine)
f.write(codeLine+"\n")
codeLine = "static void qSumC(const lbBase_t* dist, lbBase_t* ret);"
print(codeLine)
f.write(codeLine+"\n")
print("")
f.write("\n")
codeLine = "// Two phase"
print(codeLine)
f.write(codeLine+"\n")
codeLine = "static void gradPush(const lbBase_t& scalarVal, const int* neighList, VectorField<D{0:d}Q{1:d}>& grad);".format(nD, nQ)
print(codeLine)
f.write(codeLine+"\n"+"\n")
print("")
codeLine = "};"
print(codeLine)
f.write(codeLine+"\n"+"\n")
f.write("\n"+"\n")

write_dot(latticeName, nD)
write_cDot(latticeName, nD)
write_cDotAll(latticeName, nD, nQ, cBasis)
write_grad(latticeName, nD, nQ, cBasis,cLength)
write_qSum(latticeName)
write_qSumC(latticeName, nD, nQ, cBasis)
write_gradPush(latticeName, nD, nQ, cBasis,cLength)

f.write("\n")
codeLine = "#endif // LBD{0:d}Q{1:d}_H".format(nD,nQ)
print(codeLine)
f.write(codeLine+"\n")
f.write("\n")


#
f.close()


#------------------------------------------------------------------------
#----------------------------LBdXqY.cpp----------------------------------
#------------------------------------------------------------------------


f=open("LBd{0:d}q{1:d}.cpp".format(nD, nQ),"w+")

codeLine = '#include "LBd{0:d}q{1:d}.h"'.format(nD, nQ)
print(codeLine)
f.write(codeLine+"\n")
f.write("\n")

codeLine = 'constexpr lbBase_t D{0:d}Q{1:d}::w[];'.format(nD, nQ)
print(codeLine)
f.write(codeLine+"\n")
codeLine = 'constexpr int D{0:d}Q{1:d}::cDMajor_[];'.format(nD, nQ)
print(codeLine)
f.write(codeLine+"\n")
codeLine = 'constexpr lbBase_t D{0:d}Q{1:d}::cNorm[];'.format(nD, nQ)
print(codeLine)
f.write(codeLine+"\n")
codeLine = 'constexpr lbBase_t D{0:d}Q{1:d}::B[];'.format(nD, nQ)
print(codeLine)
f.write(codeLine+"\n")
f.write("\n")

f.close()




