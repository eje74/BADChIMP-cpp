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
    cl = "{"
    print(cl)
    cl = "    return "
    for d in range(nd):
        cl += "leftVec[{0:d}]*rightVec[{0:d}]".format(d)
        if d != nd-1:
            cl += " + "
        else:
            cl += ";"
    print(cl)
    cl = "}"
    print(cl)

# inline lbBase_t D2Q9::cDot(const int qDir, const lbBase_t* rightVec)
# {
#     return c(qDir, 0)*rightVec[0] + c(qDir, 1)*rightVec[1];
# }
def write_cDot(dxqy, nd):
    cl = "inline lbBase_t {0:s}::cDot(const int qDir, const lbBase_t* rightVec)".format(dxqy)
    print(cl)
    cl = "{"
    print(cl)
    cl = "    return "
    for d in range(nd):
        cl += "c(qDir, {0:d})*rightVec[{0:d}]".format(d)
        if d != nd-1:
            cl += " + "
        else:
            cl += ";"
    print(cl)
    cl = "}"
    print(cl)


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
    cl = "{"
    print(cl)
    for q in range(nq):
        cl = "ret[{0:d}] = ".format(q)
        if all([x == 0 for x in cv[q]]):
            cl += " 0.0;"
        else:
            for d in range(nd):
                cl += int_x_vec(cv[q][d], "vec", d)
            cl + ";"
        print(cl)
    cl = "}"
    print(cl)




latticeName = "D2Q9"
codeLine = "struct " + latticeName + "{"
print(codeLine)

nD = 2
nQ = 9
codeLine = "static constexpr int nD = " + str(nD) + ";"
print(codeLine)
codeLine = "static constexpr int nQ = " + str(nQ) + ";"
print(codeLine)

nDirPairs = 4
nQNonZero = 8

codeLine = "static constexpr int nDirPairs_ = " + str(nDirPairs) + ";"
print(codeLine)
codeLine = "static constexpr int nQNonZero_ = " + str(nQNonZero) + ";"
print(codeLine)

c2Inv = 3.0
c4Inv = 9.0

codeLine = "static constexpr lbBase_t c2Inv = " + str(c2Inv) + ";"
print(codeLine)
codeLine = "static constexpr lbBase_t c4Inv = " + str(c4Inv) + ";"
print(codeLine)
codeLine = "static constexpr lbBase_t c2 = 1.0 / c2Inv;"
print(codeLine)
codeLine = "static constexpr lbBase_t c4 = 1.0 / c4Inv;"
print(codeLine)
codeLine = "static constexpr lbBase_t c4Inv0_5 = 0.5 * c4Inv;"
print(codeLine)

wGcd = 36
wN = [16, 4, 1]
for num, wn in enumerate(wN):
    codeLine = "static constexpr lbBase_t w{0:d} = {1:.1f}/{2:.1f};".format(num, wn, wGcd)
    print(codeLine)
    codeLine = "static constexpr lbBase_t w{0:d}c2Inv = w{0:d}*c2Inv;".format(num)
    print(codeLine)

#### vector
# // Remember do define the static arrays in the cpp file as well.

vec = []
vecDimX = [1, 1, 0, -1, -1, -1, 0, 1, 0]
vecDimY = [0, 1, 1, 1, 0, -1, -1, -1, 0]
vec.append(vecDimX)
vec.append(vecDimY)
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
print(codeLineCVec)
print(codeLineCNorm)

# // Two phase values


bGcd = 108
bN = [-16, 8, 5]

codeLine = ""
for num, bn in enumerate(bN):
    codeLine = "static constexpr lbBase_t B{0:d} = {1:.1f}/{2:.1f};".format(num, bn, bGcd)
    print(codeLine)

codeLine = "static constexpr lbBase_t B[{0:d}] = ".format(nQ) + "{"
for q in range(nQ):
    codeLine += "B{0:d}".format(cLength[q])
    if q != nQ-1:
        codeLine += ", "
    else:
        codeLine += "};"

print(codeLine)

codeLine = "// Functions"
print(codeLine)
codeLine = "inline static int c(const int qDirection, const int dimension)  {return cDMajor_[nD*qDirection + dimension];}"
print(codeLine)
codeLine = "inline static int reverseDirection(const int qDirection) {return (qDirection + nDirPairs_) % nQNonZero_;}"
print(codeLine)
codeLine = "static lbBase_t dot(const lbBase_t* leftVec, const lbBase_t* rightVec);"
print(codeLine)
codeLine = "static lbBase_t cDot(const int qDir, const lbBase_t* rightVec);"
print(codeLine)
codeLine = "static void cDotAll(const lbBase_t* vec, lbBase_t* ret);"
print(codeLine)

write_dot(latticeName, nD)
write_cDot(latticeName, nD)
print(cBasis[8] == 0)
write_cDotAll(latticeName, nD, nQ, cBasis)
#
#
# // Functions
#
# static void cDotAll(const lbBase_t* vec, lbBase_t* ret);
# static void grad(const lbBase_t* rho, lbBase_t* ret);
#
# static void qSum(const lbBase_t* dist, lbBase_t& ret);
# static void qSumC(const lbBase_t* dist, lbBase_t* ret);
#
# // Two phase
# static void gradPush(const lbBase_t& scalarVal, const int* neighList, VectorField<D2Q9>& grad);
#
# };
#
#


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
# inline void D2Q9::grad(const lbBase_t* rho, lbBase_t* ret)
# {
#     ret[0] =  w1c2Inv * (rho[0] - rho[4]);
#     ret[0] += w2c2Inv * (rho[1] - rho[3] - rho[5] + rho[7]);
#     ret[1] =  w1c2Inv * (rho[2] - rho[6]);
#     ret[1] += w2c2Inv * (rho[1] + rho[3] - rho[5] - rho[7]);
# }
#
#
#
# inline void D2Q9::qSum(const lbBase_t* dist, lbBase_t& ret)
# {
#     ret = 0.0;
#     for (int q = 0; q < nQ; ++q)
#         ret += dist[q];
# }
#
# inline void D2Q9::qSumC(const lbBase_t* dist, lbBase_t* ret)
# {
#     ret[0] = dist[0] + dist[1]            - dist[3] - dist[4]  - dist[5]           + dist[7];
#     ret[1] =           dist[1] + dist[2]  + dist[3]            - dist[5] - dist[6] - dist[7];
# }
#
#
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
#

