from ngsolve.la import InnerProduct, MultiVector
from ngsolve import Matrix, Vector

try:
    import scipy.linalg
    from scipy import random
except:
    pass

def lobpcg(mata, matm, pre, num=1, maxit=20, printrates=True, GramSchmidt=False):

    numIterations = 0
    r = mata.CreateRowVector()

    # using multivectors for better performance
    uvecs = MultiVector(r, num)
    vecs = MultiVector(r, 2 * num)

    for v in vecs[0:num]:
        v.SetRandom()
    uvecs[:] = pre * vecs[0:num]
    lams = Vector(num * [1])

    for i in range(maxit):
        numIterations += 1

        vecs[0:num] = mata * uvecs[0:num]
        vecs[0:num] -= (matm * uvecs[0:num]).Scale(lams)

        vecs[num:2 * num] = pre * vecs[0:num]
        vecs[0:num] = uvecs[0:num]

        vecs.Orthogonalize()

        asmall = InnerProduct(vecs, mata * vecs)
        msmall = InnerProduct(vecs, matm * vecs)

        ev, evec = scipy.linalg.eigh(a=asmall, b=msmall)
        lams = Vector(ev[0:num])
        if printrates:
            print(i, ":", list(lams))

        mat = Matrix(evec[:, 0:num])
        uvecs[0:num] = vecs * mat

        if (numIterations==1):
            tmp = MultiVector(r, 2 * num)
            tmp[0:2*num] = vecs
            vecs.Extend(num)
            vecs[0:2*num] = tmp[0:2 * num]

        #todo: use span{w^i,x^i,p^i} instead of span{w^i,x^i,x^{i-1}} for better stability
        vecs[2 * num:3 * num] = vecs[0:num]

    return lams, uvecs