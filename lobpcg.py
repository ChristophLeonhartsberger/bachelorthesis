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
    vecs = MultiVector(r, 3 * num)

    for v in vecs:
        r.SetRandom()
        v.data = pre * r
    lams = Vector(num * [1])

    for i in range(maxit):
        numIterations += 1

        uvecs[0:num] = mata * vecs[0:num]
        uvecs[0:num] -= (matm * vecs[0:num]).Scale(lams)

        vecs[2*num:3 * num] = pre * uvecs[0:num]

        vecs.Orthogonalize()

        asmall = InnerProduct(vecs, mata * vecs)
        msmall = InnerProduct(vecs, matm * vecs)

        ev, evec = scipy.linalg.eigh(a=asmall, b=msmall)
        lams = Vector(ev[0:num])
        if printrates:
            print(i, ":", list(lams))

        mat = Matrix(evec[:, 0:num])
        uvecs[0:num] = vecs * mat

        #use span{u^i,u^{i-1},w^i}
        vecs[num:2*num] = vecs[0:num]
        vecs[0:num] = uvecs[0:num]

    return lams, uvecs