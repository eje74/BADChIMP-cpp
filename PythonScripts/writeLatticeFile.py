# -*- coding: utf-8 -*-
# WRITE AN LATTICE CLASS USING OLD BADCHIMP LATTICE INPUT FILES
import numpy as np

def write_code_line(cl, ofs):
    print(cl)
    ofs.write(cl + "\n")

def write_code_line_end_function(cl, ofs):
    print(cl)
    ofs.write(cl + "\n"+ "\n")
    print("")


def int_x_vec(i, v_name, v_index):
    if i == 0:
        return ""
    if i == 1:
        return " +" + v_name + "[" + str(v_index) + "]"
    if i == -1:
        return " -" + v_name + "[" + str(v_index) + "]"
    return str(i) + "*" + v_name + "[" + str(v_index) + "]"


def write_dot(dxqy, nd, ofs):
    write_code_line("template <typename T1, typename T2>", ofs)
    write_code_line("inline lbBase_t {0:s}::dot(const T1 &leftVec, const T2 &rightVec)".format(dxqy), ofs)
    write_code_line("{", ofs)
    cl = "    return "
    for d in range(nd):
        cl += "leftVec[{0:d}]*rightVec[{0:d}]".format(d)
        if d != nd-1:
            cl += " + "
        else:
            cl += ";"
    write_code_line(cl, ofs)
    write_code_line_end_function("}", ofs)

def write_cDot(dxqy, nd, ofs):
    write_code_line("template<typename T>", ofs)
    write_code_line("inline T {0:s}::cDot(const int qDir, const T* rightVec)".format(dxqy), ofs)
    write_code_line("{", ofs)
    cl = "    return "
    for d in range(nd):
        cl += "c(qDir, {0:d})*rightVec[{0:d}]".format(d)
        if d != nd-1:
            cl += " + "
        else:
            cl += ";"
    write_code_line(cl, ofs)
    write_code_line_end_function("}", ofs)

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
def write_cDotAll(dxqy, nd, nq, cv, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline std::valarray<lbBase_t> {0:s}::cDotAll(const T &vec)".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("std::valarray<lbBase_t> ret(nQ);", ofs)
    for q in range(nq):
        cl = "ret[{0:d}] =".format(q)
        if all([x == 0 for x in cv[q]]):
            cl += " 0.0;"
        else:
            for d in range(nd):
                cl += int_x_vec(cv[q][d], "vec", d)
            cl += ";"
        write_code_line(cl, ofs)
    write_code_line("return ret;", ofs)
    write_code_line_end_function("}", ofs)

# inline void D2Q9::grad(const lbBase_t* rho, lbBase_t* ret)
# {
#     ret[0] =  w1c2Inv * (rho[0] - rho[4]);
#     ret[0] += w2c2Inv * (rho[1] - rho[3] - rho[5] + rho[7]);
#     ret[1] =  w1c2Inv * (rho[2] - rho[6]);
#     ret[1] += w2c2Inv * (rho[1] + rho[3] - rho[5] - rho[7]);
# }
#
def write_grad(dxqy, nd, nq, cv, cL, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline std::valarray<lbBase_t> {0:s}::grad(const T& rho)".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("std::valarray<lbBase_t> ret(nD);", ofs)
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
        write_code_line(cl, ofs)
    write_code_line("return ret;", ofs)
    write_code_line_end_function("}", ofs)
    

#inline lbBase_t D3Q19::divGrad(const T& rho)
#{
#  lbBase_t ret;
#  ret = 2*w1c2Inv * ( + rho[0] + rho[1] + rho[2] + rho[9] + rho[10]  + rho[11] )
#    + 2*w2c2Inv * ( + rho[3] + rho[4] + rho[5] + rho[6] + rho[7] + rho[8] + rho[12] + rho[13] + rho[14] + rho[15] + rho[16] + rho[17] )
#    + 2*w0c2Inv * rho[18] - 2*c2Inv*rho[18];
#
#  return ret;
#  
#}    
#
def write_divGrad(dxqy, nd, nq, cv, cL, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline lbBase_t {0:s}::divGrad(const T& rho)".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("lbBase_t ret;", ofs)
    #for d in range(nd):
    #    cl = "ret[{0:d}] =".format(d)
    #    for clength in range(max(cL) + 1):
    #        qval = [q for q in range(nq) if cv[q][d] != 0 and cL[q] == clength]
    #        if len(qval) > 0:
    #            cl += "+ 2*w{0:d}c2Inv * (".format(clength)
    #            for qv in qval:
    #                cl += int_x_vec(cv[qv][d], " rho", qv)
    #            cl += " ) "
    #    cl += ";"
    #    write_code_line(cl, ofs)
    
    cl = "ret =".format(d)
    
    for clength in range(max(cL) + 1):
        qval = [q for q in range(nq) if cL[q] == clength]
        if len(qval) > 0:
            cl += "+ 2*"
            if clength==0:
                cl += "("
            cl += " w{0:d}c2Inv".format(clength)
            if clength==0:
                cl += " - c2Inv )"
            cl +=" * ("
            for qv in qval:
                cl += int_x_vec(1, " rho", qv)
            cl += " ) "
        #else:
        #    cl += "+ 2*( w{0:d}c2Inv - c2Inv ) * (".format(clength)    
    cl += ";"
    write_code_line(cl, ofs)        

    write_code_line("return ret;", ofs)
    write_code_line_end_function("}", ofs)
        

# inline void D2Q9::qSum(const lbBase_t* dist, lbBase_t& ret)
# {
#     ret = 0.0;
#     for (int q = 0; q < nQ; ++q)
#         ret += dist[q];
# }

def write_qSum(dxqy, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline lbBase_t {0:s}::qSum(const T &dist)".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("lbBase_t ret = 0.0;", ofs)
    write_code_line("for (int q = 0; q < nQ; ++q)", ofs)
    write_code_line("ret += dist[q];", ofs)
    write_code_line("return ret;", ofs)
    write_code_line_end_function("}", ofs)

# inline void D2Q9::qSumC(const lbBase_t* dist, lbBase_t* ret)
# {
#     ret[0] = dist[0] + dist[1]            - dist[3] - dist[4]  - dist[5]           + dist[7];
#     ret[1] =           dist[1] + dist[2]  + dist[3]            - dist[5] - dist[6] - dist[7];
# }

def write_qSumC(dxqy, nd, nq, cv, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline std::valarray<lbBase_t> {0:s}::qSumC(const T &dist)".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("std::valarray<lbBase_t> ret(nD);", ofs)
    for d in range(nd):
        cl = "ret[{0:d}] =".format(d)
        for q in range(nq):
            if any([x != 0 for x in cv[q]]):
                if(cv[q][d]!=0):
                    cl += int_x_vec(cv[q][d]," dist",q)
        write_code_line(cl + ";", ofs)
    write_code_line("return ret;", ofs)
    write_code_line_end_function("}", ofs)

def write_qSumCC(dxqy, nd, nq, cv, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline std::valarray<lbBase_t> {0:s}::qSumCCLowTri(const T &dist)".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("std::valarray<lbBase_t> ret(nD*(nD+1)/2);", ofs)
    it=0;
    for di in range(nd):
        #for dj in np.linspace(di,nd,nd-di,endpoint=False, dtype=int):
        for dj in range(di+1):    
            cl = "ret[{0:d}] =".format(it)
            for q in range(nq):
                if any([x != 0 for x in cv[q]]):
                    if(cv[q][di]!=0 and cv[q][dj]!=0):
                        cl += int_x_vec(cv[q][di]*cv[q][dj]," dist",q)
            write_code_line(cl + ";", ofs)
            it=it+1;
            
    write_code_line("return ret;", ofs)
    
    write_code_line_end_function("}", ofs)
    
def write_traceLowTri(dxqy, nd, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline lbBase_t {0:s}::traceLowTri(const T &lowTri)".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("lbBase_t ret;", ofs)
    cl = "return ret ="
    it=0;
    for di in range(nd):
        for dj in range(di+1):   
            if dj == di:
                cl += "+ lowTri[{0:d}]".format(it)
            it=it+1    
    write_code_line(cl + ";", ofs)
    write_code_line_end_function("}", ofs)
    
def write_traceOfMatrix(dxqy, nd, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline lbBase_t {0:s}::traceOfMatrix(const T &mat)".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("lbBase_t ret;", ofs)
    cl = "return ret ="
    it=0;
    for di in range(nd):
        for dj in range(nd):   
            if dj == di:
                cl += "+ mat[{0:d}]".format(it)
            it=it+1    
    write_code_line(cl + ";", ofs)
    write_code_line_end_function("}", ofs)    
    
def write_deltaMatrix(dxqy, nd, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline std::valarray<lbBase_t> {0:s}::deltaMatrix()".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("std::valarray<lbBase_t> ret(nD*nD);", ofs)
    cl = "return ret ="
    it=0;
    for di in range(nd):
        for dj in range(nd):   
            if dj == di:
                cl = "ret[{0:d}] =".format(it)+ " 1"
            else:
                cl = "ret[{0:d}] =".format(it)+ " 0"
                
            it=it+1     
            write_code_line(cl + ";", ofs)
    write_code_line_end_function("}", ofs)    
    
    
            
def write_contractionLowTri(dxqy, nd, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline lbBase_t {0:s}::contractionLowTri(const T &lowTri1, const T &lowTri2)".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("lbBase_t ret;", ofs)
    cl = "return ret ="
    it=0;
    for di in range(nd):
        for dj in range(di+1):   
            if dj == di:
                cl += "+ lowTri1[{0:d}]*lowTri2[{0:d}]".format(it)
            else:
                cl += "+ 2*lowTri1[{0:d}]*lowTri2[{0:d}]".format(it)
            it=it+1    
    write_code_line(cl + ";", ofs)
    write_code_line_end_function("}", ofs)   
    
def write_contractionRank2(dxqy, nd, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline lbBase_t {0:s}::contractionRank2(const T &mat1, const T &mat2)".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("lbBase_t ret;", ofs)
    cl = "return ret ="
    it=0;
    for di in range(nd):
        for dj in range(nd):   
            cl += "+ mat1[{0:d}]*mat2[{0:d}]".format(it)
            it=it+1    
    write_code_line(cl + ";", ofs)
    write_code_line_end_function("}", ofs)       

def write_matrixMultiplication(dxqy, nd, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline std::valarray<lbBase_t> {0:s}::matrixMultiplication(const T &mat1, const T &mat2)".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("std::valarray<lbBase_t> ret(nD*nD);", ofs)
    
    nRows=nd
    nColumns=nd
        
    for i in range(nRows):
        for j in range(nColumns):
            cl = "ret[{0:d}] =".format(j + i*nRows)
            for k in range(nColumns):        
                cl += " + mat1[{0:d}]".format(k + i*nRows)+"*mat2[{0:d}]".format(j + k*nRows)
            write_code_line(cl + ";", ofs)
            
    write_code_line("return ret;", ofs)
    write_code_line_end_function("}", ofs)             
            
def write_contractionLowTriVec(dxqy, nd, ofs):
    write_code_line("template <typename T>", ofs)
    write_code_line("inline std::valarray<lbBase_t> {0:s}::contractionLowTriVec(const T &lowTri, const T &vec)".format(dxqy), ofs)
    write_code_line("{", ofs)
    write_code_line("std::valarray<lbBase_t> ret(nD);", ofs)
    
    strArr=np.empty(nd*nd, dtype='object')
    it=0;
    for i in range(nd):
        for j in range(i+1):    
            strArr[i*nd + j] = "lowTri[{0:d}]".format(it)
            it += 1
    for i in range(nd):
        for j in np.linspace(i+1,nd, nd-(i+1),endpoint=False,dtype=int):
            strArr[i*nd + j] = strArr[j*nd + i]
    
    
    it=0
    for i in range(nd):
        cl = "ret[{0:d}] =".format(i)
        for j in range(nd):
            cl += " + "+strArr[it]+"*vec[{0:d}]".format(j)
            it += 1
        write_code_line(cl + ";", ofs)
    write_code_line("return ret;", ofs)
    write_code_line_end_function("}", ofs)         
    
    #for d in range(nd):
    #    cl = "ret[{0:d}] = 0;".format(d)    
    #for d in range(nd):
    #    cl = "ret[{0:d}] =".format(d)
    #    for q in range(nq):
    #        if any([x != 0 for x in cv[q]]):
    #            if(cv[q][d]!=0):
    #                cl += int_x_vec(cv[q][d]," dist",q)
    #    write_code_line(cl + ";", ofs)
    #write_code_line("return ret;", ofs)
    #write_code_line_end_function("}", ofs)    

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

def write_gradPush(dxqy, nd, nq, cv, cL, ofs):
    write_code_line("inline void {0:s}::gradPush(const lbBase_t &scalarVal, const int *neigList, VectorField<{0:s}> &grad)".format(dxqy), ofs)
    write_code_line("{", ofs)
    for clength in range(max(cL) + 1):
        if clength > 0:
            write_code_line("const lbBase_t valTmp{0:d}  = scalarVal * c2Inv * w{0:d};".format(clength), ofs)
    write_code_line("", ofs)

    for q in range(nq-1):
        cl = ""
        if q == 0:
            cl += "int "
        cl += "nodeNeigNo = neigList[{0:d}];".format(q)
        write_code_line(cl, ofs)
        for d in range(nd):
            if cv[q][d]!=0:
                cl = "grad(0,{0:d},nodeNeigNo) ".format(d)
                if cv[q][d]>0:
                    cl += "-= "
                else:
                    cl += "+= "
                cl += "valTmp{0:d};".format(cL[q])
                write_code_line(cl, ofs)
        write_code_line("", ofs)
    write_code_line_end_function("}", ofs)


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


# WRITE_FILE_HEADER
f=open("LBd{0:d}q{1:d}.h".format(nD, nQ),"w+")

write_code_line("#ifndef LBD{0:d}Q{1:d}_H".format(nD,nQ), f)
write_code_line("#define LBD{0:d}Q{1:d}_H".format(nD,nQ) + "\n", f)

# -- write include
write_code_line('#include "LBglobal.h"', f)
write_code_line('#include "LBfield.h"', f)
write_code_line('#include <vector>' + "\n", f)
write_code_line('// See "LBlatticetypes.h" for description of the structure' + "\n", f)

# WRITE STRUCT
latticeName = "D{0:d}Q{1:d}".format(nD,nQ)
write_code_line("struct " + latticeName + "{" + "\n", f)
write_code_line("static constexpr int nD = " + str(nD) + ";", f)
write_code_line("static constexpr int nQ = " + str(nQ) + ";", f)
write_code_line("static constexpr int nDirPairs_ = " + str(nDirPairs) + ";", f)
write_code_line("static constexpr int nQNonZero_ = " + str(nQNonZero) + ";" + "\n", f)
write_code_line("static constexpr lbBase_t c2Inv = " + str(c2Inv) + ";", f)
write_code_line("static constexpr lbBase_t c4Inv = " + str(c4Inv) + ";", f)
write_code_line("static constexpr lbBase_t c2 = 1.0 / c2Inv;", f)
write_code_line("static constexpr lbBase_t c4 = 1.0 / c4Inv;", f)
write_code_line("static constexpr lbBase_t c4Inv0_5 = 0.5 * c4Inv;" + "\n", f)

for num, wn in enumerate(wN):
    write_code_line("static constexpr lbBase_t w{0:d} = {1:.1f}/{2:.1f};".format(num, wn, wGcd), f)
    write_code_line("static constexpr lbBase_t w{0:d}c2Inv = w{0:d}*c2Inv;".format(num), f)
write_code_line("", f)

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
write_code_line(codeLineW, f)
write_code_line(codeLineCVec, f)
write_code_line(codeLineCNorm, f)

codeLineRevDir = "static constexpr int reverseDirection_[{0:d}] = ".format(nQ) + "{"

for q in range(nQNonZero):
    codeLineRevDir += "{0:d}, ".format((q + nDirPairs)%nQNonZero)
codeLineRevDir += "{0:d}".format(nQNonZero)
codeLineRevDir += "};"
write_code_line(codeLineRevDir, f)

# // Two phase values
for num, bn in enumerate(bN):
    codeLine = "static constexpr lbBase_t B{0:d} = {1:.1f}/{2:.1f};".format(num, bn, bGcd)
    write_code_line(codeLine, f)

codeLine = "static constexpr lbBase_t B[{0:d}] = ".format(nQ) + "{"
for q in range(nQ):
    codeLine += "B{0:d}".format(cLength[q])
    if q != nQ-1:
        codeLine += ", "
    else:
        codeLine += "};"
write_code_line(codeLine+"\n", f)
codeLine = "static constexpr lbBase_t UnitMatrixLowTri[{0:d}] = ".format(nD*(nD+1)/2) + "{"
it=0
for di in range(nD):
    for dj in range(di+1):
        if dj == di:
            codeLine += "1"
        else:
            codeLine += "0"
        it=it+1    
        if it != (nD*(nD+1)/2):
           codeLine += ", "        
        else:
            codeLine += "};"
write_code_line(codeLine+"\n", f)


write_code_line("// Functions" + "\n", f)
write_code_line("inline static int c(const int qDirection, const int dimension)  {return cDMajor_[nD*qDirection + dimension];}", f)
#write_code_line("inline static int reverseDirection(const int qDirection) {return (qDirection + nDirPairs_) % nQNonZero_;}" + "\n", f)
write_code_line("inline static int reverseDirection(const int qDirection) {return reverseDirection_[qDirection];}" + "\n", f)

write_code_line("template <typename T1, typename T2>", f)
write_code_line("inline static lbBase_t dot(const T1 &leftVec, const T2 &rightVec);", f)

write_code_line("template<typename T>", f)
write_code_line("inline static T cDot(const int qDir, const T* rightVec);", f)

write_code_line("template <typename T>", f)
write_code_line("inline static std::valarray<lbBase_t> cDotAll(const T &vec);", f)
write_code_line("template <typename T>", f)
write_code_line("inline static std::valarray<lbBase_t> grad(const T &rho);" + "\n", f)
write_code_line("template <typename T>", f)
write_code_line("inline static lbBase_t divGrad(const T &rho);" + "\n", f)

write_code_line("template <typename T>", f)
write_code_line("inline static lbBase_t qSum(const T &dist);", f)
write_code_line("template <typename T>", f)
write_code_line("inline static std::valarray<lbBase_t> qSumC(const T &dist);" + "\n", f)
write_code_line("template <typename T>", f)
write_code_line("inline static std::valarray<lbBase_t> qSumCCLowTri(const T &dist);" + "\n", f)
write_code_line("template <typename T>", f)
write_code_line("inline static lbBase_t traceLowTri(const T &lowTri);" + "\n", f)
write_code_line("template <typename T>", f)
write_code_line("inline static lbBase_t traceOfMatrix(const T &mat);" + "\n", f)
write_code_line("template <typename T>", f)
write_code_line("inline static std::valarray<lbBase_t> deltaMatrix();" + "\n", f)
write_code_line("template <typename T>", f)
write_code_line("inline static lbBase_t contractionLowTri(const T &lowTri1, const T &lowTri2);" + "\n", f)
write_code_line("template <typename T>", f)
write_code_line("inline static lbBase_t contractionRank2(const T &mat1, const T &mat2);" + "\n", f)
write_code_line("template <typename T>", f)
write_code_line("inline static std::valarray<lbBase_t> matrixMultiplication(const T &mat1, const T &mat2);" + "\n", f)
write_code_line("template <typename T>", f)
write_code_line("inline static std::valarray<lbBase_t> contractionLowTriVec(const T &lowTri, const T &vec);" + "\n", f)

write_code_line("// Two phase", f)
write_code_line("static void gradPush(const lbBase_t& scalarVal, const int* neighList, VectorField<D{0:d}Q{1:d}>& grad);".format(nD, nQ) + "\n", f)
write_code_line("};" + "\n"+"\n", f)

# WRITE FUNCTION DEFINITIONS
write_dot(latticeName, nD, f)
write_cDot(latticeName, nD, f)
write_cDotAll(latticeName, nD, nQ, cBasis, f)
write_grad(latticeName, nD, nQ, cBasis,cLength, f)
write_divGrad(latticeName, nD, nQ, cBasis,cLength, f)
write_qSum(latticeName, f)
write_qSumC(latticeName, nD, nQ, cBasis, f)
write_qSumCC(latticeName, nD, nQ, cBasis, f)
write_traceLowTri(latticeName, nD, f)
write_traceOfMatrix(latticeName, nD, f)
write_deltaMatrix(latticeName, nD, f)
write_contractionLowTri(latticeName, nD, f)
write_contractionRank2(latticeName, nD, f)
write_matrixMultiplication(latticeName, nD, f)
write_contractionLowTriVec(latticeName, nD, f)
write_gradPush(latticeName, nD, nQ, cBasis, cLength, f)

write_code_line("", f)
write_code_line("#endif // LBD{0:d}Q{1:d}_H".format(nD,nQ), f)


#
f.close()


#------------------------------------------------------------------------
#----------------------------LBdXqY.cpp----------------------------------
#------------------------------------------------------------------------


f=open("LBd{0:d}q{1:d}.cpp".format(nD, nQ),"w+")

write_code_line('#include "LBd{0:d}q{1:d}.h"'.format(nD, nQ) + "\n", f)

write_code_line('constexpr lbBase_t D{0:d}Q{1:d}::w[];'.format(nD, nQ), f)
write_code_line('constexpr int D{0:d}Q{1:d}::cDMajor_[];'.format(nD, nQ), f)
write_code_line('constexpr lbBase_t D{0:d}Q{1:d}::cNorm[];'.format(nD, nQ), f)
write_code_line('constexpr int D{0:d}Q{1:d}::reverseDirection_[];'.format(nD, nQ), f)
write_code_line('constexpr lbBase_t D{0:d}Q{1:d}::B[];'.format(nD, nQ) + "\n", f)

f.close()
