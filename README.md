<h1 align="center">
  Response-time analysis for preemptive jobs using SAG
</h1>
<h4 align="center">Reachability-Based Response-Time Analysis of Preemptive Tasks Under Global Scheduling using Schedule Abstraction Graph (SAG)</h4>

This repository contains the implementations of response-time analysis for **sets of preemptive jobs** scheduled *
*globally on a system with identical cores**.
The analysis is described in the following paper:

- P. Gohari, J. Voeten, and M. Nasri “Reachability-Based Response-Time Analysis of Preemptive Tasks Under Global
  Scheduling”, *Euromicro Conference on Real-Time Systems (ECRTS)*, 2024.

## 📦 Dependencies

- A modern C++ compiler supporting the **C++14 standard**. Recent versions of `clang` and `g++` on Linux is known to
  work.
- The [CMake](https://cmake.org) build system. For installation using `apt` (Ubuntu, Debian...):

```bash
sudo apt-get -y install cmake 
``` 

- A POSIX OS. Linux and macOS are both known to work.

- The [Intel Thread Building Blocks (TBB)](https://www.threadingbuildingblocks.org) library and parallel runtime.

- The [jemalloc](http://jemalloc.net) scalable memory allocator. Alternatively, the TBB allocator can be used instead;
  see build options below.

## 📋 Build Instructions

These instructions assume a Linux host.

To compile the tool, first generate an appropriate `Makefile` with `cmake` and then use it to actually build the source
tree.

```bash
 # (1) enter the build directory
 cd build
 # (2) generate the Makefile
 cmake ..
 # (3) build everything
 make -j
```

The last step yields the main binary `rt_analysis` in the `build/` directory.

### Build Options

The build can be configured in a number of ways by passing options via the `-D` flag to `cmake` in step (2).

To enable debug builds, pass the `DEBUG` option to `cmake` .

    cmake -DDEBUG=yes ..

To enable the collection of schedule graphs (the `-g` option in `rt_analysis`), set the option `COLLECT_SCHEDULE_GRAPHS`
to `yes`.

    cmake -DCOLLECT_SCHEDULE_GRAPHS=yes  ..

Note that enabling `COLLECT_SCHEDULE_GRAPHS` turns off parallel analysis, i.e., the analysis becomes single-threaded, so
don't turn it on by default. It is primarily a debugging aid.

By default, `rt_analysis` uses the default `libc` allocator (which is a tremendous scalability bottleneck).
To use the `jemalloc` or `Intel TBB` allocator instead, pass the `USE_JE_MALLOC` or `USE_TBB_MALLOC` option to `cmake`.

## 📝 Input Format

The tool operates on classical SAG CSV files with a fixed column order. For more details, see the examples in
the `examples/` folder.

Job set input CSV files describe a set of jobs, where each row specifies one job. The following columns are required.

1. **Task ID** — an arbitrary numeric ID to identify the task to which a job belongs
2. **Job ID** — a unique numeric ID that identifies the job
3. **Release min** — the earliest-possible release time of the job (equivalently, this is the arrival time of the job)
4. **Release max** — the latest-possible release time of the job (equivalently, this is the arrival time plus maximum
   jitter of the job)
5. **Cost min** — the best-case execution time of the job (can be zero)
6. **Cost max** — the worst-case execution time of the job
7. **Deadline** — the absolute deadline of the job
8. **Priority** — the priority of the job (EDF: set it equal to the deadline)

All numeric parameters can be 64-bit integers (preferred) or floating point values (slower, not recommended).

## ⚙️ Usage

To run the tool on a given set, pass the filename as an argument and provide the number of (identical) processors (cores) via the `-m` option.  
For example, to analyze the job set in [examples/fig1-pr.csv](examples/fig1-pr.csv) for a system with two cores, run the
following command:

```bash
./build/rt_analysis examples/fig1-pr.csv -m 2
examples/fig1-pr.csv,  1,  5,  24,  30,  7,  0.000044,  6.023438,  0,  2
```

See the builtin help (`rt_analysis -h`) for further options.

## 📄 Output Format

The output is provided in CSV format and consists of the following columns:

1. The input file name.
2. The schedulability result:
    - 1 if the job is *is* schedulable (i.e., the tool could prove the absence of deadline misses),
    - 0 if it is *not*, or if the analysis timed out, if it reached the depth limit, or if the analysis cannot prove the
      absence of deadline misses (while the RTSS'17 analysis is exact, the ECRTS'19 analysis is only sufficient, but not
      exact).
3. The number of jobs in the job set.
4. The number of states that were explored.
5. The number of edges that were discovered.
6. The maximum “exploration front width” of the schedule graph, which is the maximum number of unprocessed states that
   are queued for exploration (at any point in time).
7. The CPU time used in the analysis (in seconds).
8. The peak amount of memory used (as reported by `getrusage()`), divided by 1024. Due to non-portable differences
   in `getrusage()`, on Linux this reports the memory usage in megabytes, whereas on macOS it reports the memory usage
   in kilobytes.
9. A timeout indicator: 1 if the state-space exploration was aborted due to reaching the time limit (as set with
   the `-l` option); 0 otherwise.
10. The number of processors assumed during the analysis.

Pass the `--header` flag to `rt_analysis` to print out column headers.

### Obtaining Response Times

The analysis computes for each job the earliest and latest possible completion times, from which it is trivial to infer
minimum and maximum response times. To obtain this information, pass the `-r` option to `rt_analysis`.

If invoked on an input file named `foo.csv`, the completion times will be stored in a file `foo.rta.csv` and follow the
following format, where each row corresponds to one job in the input job set.

1. Task ID
2. Job ID
3. BCCT, the best-case completion time
4. WCCT, the worst-case completion time
5. BCRT, the best-case response time (relative to the minimum release time)
6. WCRT, the worst-case response time (relative to the minimum release time)

Note that the analysis by default aborts after finding the first deadline miss, in which case some of the rows may
report nonsensical default values.

## ⚠️ Notes

- The parallel version of the analysis is implemented using Intel's Thread Building Blocks (TBB) library.
  The analysis is parallelized at the level of the state-space exploration, where each thread processes a different
  state.
  However, the analysis is not tested with parallel option enabled.
  **Use it at your own risk**.

## 🌱 Contribution

With your feedback and conversation, you can assist me in developing this framework.

* Open pull request with improvements
* Discuss feedback and bugs in issues

## 📜 License

The code is released under a 3-clause BSD license.

## 💡 Credits

Schedule Abstraction Graph (SAG) was originally introduced by Nasri et al. [1][2][3].
This implementation made use of the original tool which was developed
by [Björn Brandenburg](https://people.mpi-sws.org/~bbb/). It is now being maintained
by [Geoffrey Nelissen](https://www.tue.nl/en/research/researchers/geoffrey-nelissen/).
For more information, please visit [Schedule Abstraction Framework page](https://github.com/SAG-org).

## 📚 References

* [1] M. Nasri and B.
Brandenburg, “[An Exact and Sustainable Analysis of Non-Preemptive Scheduling](https://people.mpi-sws.org/~bbb/papers/pdf/rtss17.pdf)”, *Proceedings of the 38th IEEE Real-Time Systems Symposium (RTSS 2017)*, pp. 12–23, December 2017.

* [2] M. Nasri, G. Nelissen, and B. Brandenburg, “[A Response-Time Analysis for Non-Preemptive Job Sets under Global Scheduling](http://drops.dagstuhl.de/opus/volltexte/2018/8994/pdf/LIPIcs-ECRTS-2018-9.pdf)”, *Proceedings of the 30th Euromicro Conference on Real-Time Systems (ECRTS 2018)*, pp. 9:1–9:23, July 2018.

* [3] M. Nasri, G. Nelissen, and B. Brandenburg, “[Response-Time Analysis of Limited-Preemptive Parallel DAG Tasks under Global Scheduling](http://drops.dagstuhl.de/opus/volltexte/2019/10758/pdf/LIPIcs-ECRTS-2019-21.pdf)”, *Proceedings of the 31st Euromicro Conference on Real-Time Systems (ECRTS 2019)*, pp. 21:1–21:23, July 2019.