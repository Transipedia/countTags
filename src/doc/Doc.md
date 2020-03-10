#  Read fastq faster

 * [5 solutions](Fastq_Faster.md)

# Use Parallele

 * see some solution in previous part
 * see [c++ threads](cppthread.md)

## Docs
 * [https://attractivechaos.wordpress.com/2013/10/11/parallelizing-simple-for-loops/](simple loops C)
 * [http://attractivechaos.github.io/klib/#Kthread%3A%20simple%20threading%20models](kthreads from klib)
 * [https://bisqwit.iki.fi/story/howto/openmp/](guide to open mp)
 * [https://people.sc.fsu.edu/~jburkardt/cpp_src/openmp/openmp.html](example with open mp)
 * [https://stackoverflow.com/questions/2352895/how-to-ensure-a-dynamically-allocated-array-is-private-in-openmp?noredirect=1](private array in open mp)


# Libraries to use may be

## Benchmark

 * [google benchmar](https://github.com/google/benchmark.git): easy to use (define macro)

## DB storing

 * [use of badger (GO) to store all kmer](https://github.com/dgraph-io/badger)

# Usefull code not needed for now

##  get info from enviroment variable
    #include <cstdlib>

    // manage VERSION via git and GIT_VERSION environment variable
      std::string VERSION = "unrelease";
      char * val = std::getenv("GIT_VERSION");
      if (val != NULL)
        VERSION = std::string(val);

