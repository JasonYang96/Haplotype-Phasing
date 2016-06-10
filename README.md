# Haplotype-Phasing
## CS CM124 Computational Genetics
Jason Yang 804331785

# Description of the Project
Current sequencing technologies, although very cost-effective and high-throughput, can not distinguish the haplotype a parent passed on to their children. There are many applications of haplotype phasing as studies of haplotype lineages can lead to mre discoveries about how genes and traits are passed on as well as how genetic diseases work in terms of haplotypes.

Computationally, this program will take in several varying difficulties of genotype reads and produce the minimum optimal set of haplotypes. I have coded the trivial solution as well as the improved solution to compare. The benchmarks I will be using are computational time to produce the minimal set, the number of haplotypes produced by the algorithm, and memory requirements.

The trivial solution is the Cartesian product of all possible ambiguous SNPs in the genotype. The improved solutions is essentially Clark's Method with a slight improvement. My goal was to decrease the number of unphased genotypes, so before I run another iteration of Clark's Method, I check to see if any unphased genotypes can be made by combining two known haplotypes.

What I have learned from this project is that the improved solution beats the trivial solution in every benchmark. In terms of the "improvement" made to Clark's Method, going through varying difficulties of data, the algorithm has shown that this improvement is slight at best. It depends on a number of varying factors and will only mostly be an improvement if there are considerably more unphased genotypes than known haplotypes as it will be worth it to compute the combinations of known genotypes than to run another iteration of Clark's Method.

From here, improvement begins by continuing to look at Clark's Method to see if there are any places where I can reduce the number of unphased genotypes with known haplotypes. I want to reduce this set of unphased genotypes as much as possible as each iteration of Clark's Method is computationally heavy. Two places to look would be smarter selection of known haplotypes to combine with unphased genotypes, and smarter creation of new haplotypes when phasing genotypes.

# Description of Me
My name is Jason Yang and I a just finished my Junior Year at the University of California, Los Angeles. I'm currently an undergrad studying Computer Science and mainly took this class as a computer science elective.

# Goal for end of the quarter
My goal for the end of the quarter is to be able to at least complete the medium difficulty of the project description, which is to be able to phase up to 15 SNPs as well as compute minimum parsimony haplotypes for at least 10 SNPs

#Weekly Schedule
- Week 1: looked over project
- Week 2: rewatched lecture on project 1
- Week 3: started coding trivial and improved solutions
- Week 4: signed up for a presentation and finished cleaning up project
