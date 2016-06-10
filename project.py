# uses python 2, seems to run into undefined behavior with python 3
import numpy as np
import timeit

### load file as np array
def load(file_name):
    genotype = []
    f = open(file_name, 'rb')
    for line in f:
        n = list(line.strip())
        genotype.append(n)
    f.close()
    genotype = np.array(genotype, dtype=np.int8)
    return genotype

### sort their haplotypes to compare
def check(theirs, str):
    theirs = unique_rows(theirs)
    np.savetxt(str,theirs, fmt = '%s', delimiter="")

### sorts rows of a and removes duplicate rows
def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

### produces cartesian product of arrays
def cartesian(arrays, out=None):
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in range(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

### Takes a matrix, 'a', of genotype reads
### and saves the resulting haplotype
### phases to the file 'str'
### provides all possible haplotypes
def trivial_method(a, str):
    haplotypes = []
    haplotype = []
    out = np.empty([1, a.shape[1]])
    for x in range(len(a)):
        line = a[x]
        for i in range(len(line)):
            if line[i] == 1:
                haplotype.append([0,1])
            elif line[i] == 0:
                haplotype.append([0])
            else:
                haplotype.append([1])
        haplotypes.append(haplotype)
        haplotype = []
    haplotypes = np.array(haplotypes)
    for haplotype in haplotypes:
        line = cartesian(haplotype)
        out = np.append(out, line, axis=0)
    out = np.delete(out, 0, 0)
    out = np.array(out, dtype=np.int8)
    out = unique_rows(out)
    np.savetxt(str, out, fmt = '%s', delimiter="")

### Takes a matrix, 'a', of genotype reads
### and saves the resulting haplotype
### phasing to the file 'str'
### Provides optimum minimum parsimony
def imp_method(a, str):
    haplotypes = []
    haplotype = []
    delete = []
    # first get known haplotypes
    for x in range(len(a)):
        line = a[x]
        for i in range(len(line)):
            if line[i] == 1:
                break
            elif line[i] == 0:
                haplotype.append(0)
            else:
                haplotype.append(1)
        if len(haplotype) == len(line):
            haplotypes.append(haplotype)
            delete.append(x)
        haplotype = []    
    haplotypes = np.array(haplotypes)
    haplotypes = unique_rows(haplotypes)
    a = np.delete(a, delete, 0)

    # then go through improved clark's method's iterations
    while len(a) > 0:
        for haplotype in haplotypes:
            diff = a - haplotype
            zero_rows = np.where(np.all(diff >= 0,axis=1))[0]
            two_rows = np.where(np.all(diff < 2,axis=1))[0]
            rows = np.intersect1d(zero_rows, two_rows)
            haplotypes = np.append(haplotypes, diff[rows], axis=0)
            haplotypes = unique_rows(haplotypes)
            a = np.delete(a, rows, 0)
    np.savetxt(str,haplotypes, fmt = '%s', delimiter="")

# load genotypes
veasy = load("data/very_easy_training_genotypes.txt")
easy = load("data/easy_training_genotypes.txt")
med = load("data/medium_training_genotypes.txt")
hard = load("data/hard_training_genotypes.txt")
veasy_theirs = load("data/very_easy_training_haplotypes.txt")
easy_theirs = load("data/easy_training_haplotypes.txt")
medium_theirs = load("data/medium_training_haplotypes.txt")
hard_theirs = load("data/hard_training_haplotypes.txt")
veasy_test = load("data/very_easy_test_genotypes.txt")
easy_test = load("data/easy_test_genotypes.txt")
med_test = load("data/medium_test_genotypes.txt")
hard_test = load("data/hard_test_genotypes.txt")

#load haplotypes
check(veasy_theirs, "theirs_veasy.txt")
check(easy_theirs, "theirs_easy.txt")
check(medium_theirs, "theirs_medium.txt")
check(hard_theirs, "theirs_hard.txt")

# start test data
start_time = timeit.default_timer()
imp_method(veasy_test, "veasy.txt")
print("very easy test data took " + str((timeit.default_timer() - start_time) * 1000) + " ms")
start_time = timeit.default_timer()
imp_method(easy_test, "easy.txt")
print("easy test data took " + str((timeit.default_timer() - start_time) * 1000) + " ms")
start_time = timeit.default_timer()
imp_method(med_test, "med.txt")
print("medium test data took " + str((timeit.default_timer() - start_time) * 1000) + " ms")
start_time = timeit.default_timer()
imp_method(hard_test, "hard.txt")
print("hard test data took " + str((timeit.default_timer() - start_time) * 1000) + " ms")

# start improved
start_time = timeit.default_timer()
imp_method(veasy, "mine_veasy.txt")
print("very easy data took " + str((timeit.default_timer() - start_time) * 1000) + " ms")
start_time = timeit.default_timer()
imp_method(easy, "mine_easy.txt")
print("easy data took " + str((timeit.default_timer() - start_time) * 1000) + " ms")
start_time = timeit.default_timer()
imp_method(med, "mine_medium.txt")
print("medium data took " + str((timeit.default_timer() - start_time) * 1000) + " ms")
start_time = timeit.default_timer()
imp_method(hard, "mine_hard.txt")
print("hard data took " + str((timeit.default_timer() - start_time) * 1000) + " ms")

# start trivial
start_time = timeit.default_timer()
trivial_method(veasy, "trivial_veasy.txt")
print("trivial very easy data took " + str((timeit.default_timer() - start_time) * 1000) + " ms")
start_time = timeit.default_timer()
trivial_method(easy, "trivial_easy.txt")
print("trivial easy data took " + str((timeit.default_timer() - start_time) * 1000) + " ms")
# I don't have enouogh computing power to time these trivial problems
#start_time = timeit.default_timer()
#trivial_method(med, "trivial_medium.txt")
#print("trivial medium data took " + str((timeit.default_timer() - start_time) * 1000) + " ms")
# start_time = timeit.default_timer()
# trivial_method(hard, "trivial_hard.txt")
# print("trivial hard data took " + str((timeit.default_timer() - start_time) * 1000) + " ms")

