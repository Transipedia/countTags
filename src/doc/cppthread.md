# C++ 11 pool thread class
 * code at [ThreadPool](https://github.com/progschj/ThreadPool)
 * from labex.io (account labex236644; anthony.boureux@yopmail.fr/medline18), C++ thread course:

## Implemetation
        # Thread pool library with C&#43;&#43; within 100 lines: Implementation

## 1. Introduction

### Project Description

In server-side development, we usually facing the tasks of scheduling and managing the execution performance of process. In this project, we will implement a simple thread pool library by using C&#43;&#43;11.

## Things to Learn

- C&#43;&#43;11 Library Features
  - std::thread
  - std::mutex, std::unique_lock
  - std::condition_variable
  - std::future, std::packaged_task
  - std::function, std::bind
  - std::shared_ptr, std::make_shared
  - std::move, std::forward

- C&#43;&#43;11 Language Features
  - Lambda Expression
  - Trail Return Type

- Thread Pool Model
- Test Driven

## To Readers

This project is an advanced project even the code has only lower than 100 lines. Please be much more patient for it.

## 2. Working Principles

### Thread Pool Model


![image desc](/upload/M/V/Q/RFHzu4iXVYQv.png)


In the previous experiment, we learned about the thread pool model, and also reviewed many of the new features of C&#43;&#43;11 related to concurrent programming. In this section we will start to implement this thread pool within 100 lines.

### Test Driven

Thread pools are usually written as modules within a system (that is, higher-level ones can be called libraries). In this case, we must inevitably ensure the robustness of the code (at least to ensure the correctness of the library). To do this, the best way is to write the test code first, and then drive the development of the entire library based on the test code.

We want to design a `ThreadPool` class to manage the thread pool, so we can try to write such a main function to test:

```cpp
#include &lt;iostream&gt; // std::cout, std::endl

#include &lt;vector&gt;   // std::vector
#include &lt;string&gt;   // std::string
#include &lt;future&gt;   // std::future
#include &lt;thread&gt;   // std::this_thread::sleep_for
#include &lt;chrono&gt;   // std::chrono::seconds

#include &#34;ThreadPool.hpp&#34;

int main()
{
    // create a thread pool with max. 4 concurrency threads
    ThreadPool pool(4);
    // create execution results list
    std::vector&lt; std::future&lt;std::string&gt; &gt; results;

    // start eight thread task
    for(int i = 0; i &lt; 8; &#43;&#43;i) {
        // add all task to result list
        results.emplace_back(
            // ass print task to thread pool
            pool.enqueue([i] {
                std::cout &lt;&lt; &#34;hello &#34; &lt;&lt; i &lt;&lt; std::endl;
                // wait a sec when the previous line is out
                std::this_thread::sleep_for(std::chrono::seconds(1));
                // keep output and return the status of execution
                std::cout &lt;&lt; &#34;world &#34; &lt;&lt; i &lt;&lt; std::endl;
                return std::string(&#34;---thread &#34;) &#43; std::to_string(i) &#43; std::string(&#34; finished.---&#34;);
            })
        );
    }

    // outputs
    for(auto &amp;&amp; result: results)
        std::cout &lt;&lt; result.get() &lt;&lt; &#39; &#39;;
    std::cout &lt;&lt; std::endl;
    
    return 0;
}
```

### Design

We first determine the memeber and method of ThreadPool class:

```cpp

class ThreadPool {
public:

    // initialize the number of concurrency threads
    ThreadPool(size_t threads);

    // enqueue new thread task
    template&lt;class F, class... Args&gt;
    auto enqueue(F&amp;&amp; f, Args&amp;&amp;... args)
        -&gt; std::future&lt;typename std::result_of&lt;F(Args...)&gt;::type&gt;;

    // destroy thread pool and all created threads
    ~ThreadPool();
private:

    // thread list, stores all threads
    std::vector&lt; std::thread &gt; workers;

    // queue task, the type of queue elements are functions with void return type
    std::queue&lt; std::function&lt;void()&gt; &gt; tasks;

    // for synchonization
    std::mutex queue_mutex;
    // std::condition_variable is a new feature from c&#43;&#43;11,
    // it&#39;s a synchronization primitives. it can be used 
    // to block a thread or threads at the same time until
    // all of them modified condition_variable.
    std::condition_variable condition;

    // for stop
    bool stop;
};
```

In the above code, the most difficult definition is:

```cpp
template&lt;typename F, typename... Args&gt;
    auto enqueue(F&amp;&amp; f, Args&amp;&amp;... args) 
        -&gt; std::future&lt;typename std::result_of&lt;F(Args...)&gt;::type&gt;;
```

In this definition, a variable-parameter template is designed. We use `typename...` to indicate that the following parameter represents a list of zero or more types and use `Args` to represent it. The `Args` is a package of template parameter. In the later parameter list, `args` is a function parameter package that used to represent zero or more parameters.

The enqueue function is designed so that a new thread is passed in. In the first argument, use `F&amp;&amp;` to make the rvalue reference. Instead of copying the behavior, create an implementation again. Its return type, we design To obtain the result to be executed after the execution of the function body is completed, it is necessary to infer what the return type is to be able to write the code of `std::future&lt;typename return_type&gt;`.

In fact, the return type of the execution function can be inferred by `std::result_of`. We only need to use the type `F` of the function to be tuned, and the argument of this function `Args` as `std::result_of&lt;F(Args...)&gt;::type`.

Finally, this definition uses the C&#43;&#43;11 feature of the tail return type, but it is not because of the inability to write the return type, but because the return type is too long, using tails can make our code better than the following Writing is more clear and good:

```cpp
// not recommended, top-heavy
template&lt;class F, class... Args&gt;
std::future&lt;typename std::result_of&lt;F(Args...)&gt;::type&gt; 
enqueue(F&amp;&amp; f, Args&amp;&amp;... args)
```

## 3. Implementation

In the familiarity with language features, concurrent programming basic knowledge, and a deep understanding of the thread pool model, coding is not complicated anymore.

### Construction and deconstruction

We first implement the thread pool constructor. In the process of writing a thread pool constructor, pay attention to the following points:

1. The constructed thread pool should perform the specified number of concurrently executable threads;
2. The task&#39;s execution and completion phase should be in the critical section, so that we can ensure concurrently starting multiple tasks that need to be performed will not occur concurrently;

Here is an implementation:

```cpp
// constructor initialize a fixed size of worker
inline ThreadPool::ThreadPool(size_t threads): stop(false)
{

    // initialize worker
    for(size_t i = 0;i&lt;threads;&#43;&#43;i)
        // std::vector::emplace_back :
        //    append to the end of vector container
        //    this element will be constructed at the end of container, without copy and move behavior
        workers.emplace_back(
            // the lambda express capture this, i.e. the instance of thread pool
            [this]
            {
                // avoid fake awake
                for(;;)
                {

                    // define function task container, return type is void
                    std::function&lt;void()&gt; task;

                    // critical section
                    {
                        // get mutex
                        std::unique_lock&lt;std::mutex&gt; lock(this-&gt;queue_mutex);

                        // block current thread
                        this-&gt;condition.wait(lock,
                            [this]{ return this-&gt;stop || !this-&gt;tasks.empty(); });

                        // return if queue empty and task finished
                        if(this-&gt;stop &amp;&amp; this-&gt;tasks.empty())
                            return;

                        // otherwise execute the first element of queue
                        task = std::move(this-&gt;tasks.front());
                        this-&gt;tasks.pop();
                    }

                    // execution
                    task();
                }
            }
        );
}
```

The destruction of the thread pool corresponds to what instance was created during the construction.
Before destroying the thread pool, we do not know whether the execution of the worker thread in the current thread pool is completed.
Therefore, a critical section must be created to mark the thread pool as stopped. The new thread joins, and finally waits for the execution of all execution threads to finish and complete the destruction. The detailed implementation is as follows:

```cpp
// destroy everything
inline ThreadPool::~ThreadPool()
{
    // critical section
    {
        std::unique_lock&lt;std::mutex&gt; lock(queue_mutex);
        stop = true;
    }
    
    // wake up all threads
    condition.notify_all();
    
    // let all processes into synchronous execution, use c&#43;&#43;11 new for-loop: for(value:values)
    for(std::thread &amp;worker: workers)
        worker.join();
}
```

### Enqueue for new tasks

Adding a new task&#39;s implementation logic to the thread pool requires the following major considerations:

1. Support variable length template parameters for multiple enqueue task parameters
2. In order to schedule the task to be performed, the tasks to be executed need to be packaged, which means that the type of the task function needs to be packaged and constructed.
3. The critical section can be created within a scope. The best practice is to use the RAII form

```cpp
// Enqueue a new thread
// use variadic templates and tail return type
template&lt;class F, class... Args&gt;
auto ThreadPool::enqueue(F&amp;&amp; f, Args&amp;&amp;... args)
    -&gt; std::future&lt;typename std::result_of&lt;F(Args...)&gt;::type&gt;
{
    // deduce return type
    using return_type = typename std::result_of&lt;F(Args...)&gt;::type;

    // fetch task
    auto task = std::make_shared&lt; std::packaged_task&lt;return_type()&gt; &gt;(
        std::bind(std::forward&lt;F&gt;(f), std::forward&lt;Args&gt;(args)...)
    );

    std::future&lt;return_type&gt; res = task-&gt;get_future();

    // critical section
    {
        std::unique_lock&lt;std::mutex&gt; lock(queue_mutex);

        // avoid add new thread if theadpool is destroied
        if(stop)
            throw std::runtime_error(&#34;enqueue on stopped ThreadPool&#34;);

        // add thread to queue
        tasks.emplace([task]{ (*task)(); });
    }

    // notify a wait thread
    condition.notify_one();
    return res;
}
```

### Result

Compile the code and execute it:

```bash
g&#43;&#43; main.cpp -std=c&#43;&#43;11 -pthread
```

One possible result would be like this (obviously, your result is almost impossible to be exactly the same as the result here):

![image desc](/upload/D/P/A/IdBM0cnUi3fj.png)


## Conclusions

In this project, we used massive C&#43;&#43;11 feature and implemented a thread pool within 100 lines:

```cpp
// 
// ThreadPool.hpp
// ThreadPool
// 
// Original Author: Jakob Progsch, Václav Zeman
// Modified By:     changkun at https://labex.io
// Original Link:   https://github.com/progschj/ThreadPool
//

#ifndef ThreadPool_hpp
#define ThreadPool_hpp
#include &lt;vector&gt;               // std::vector
#include &lt;queue&gt;                // std::queue
#include &lt;memory&gt;               // std::make_shared
#include &lt;stdexcept&gt;            // std::runtime_error
#include &lt;thread&gt;               // std::thread
#include &lt;mutex&gt;                // std::mutex,        std::unique_lock
#include &lt;condition_variable&gt;   // std::condition_variable
#include &lt;future&gt;               // std::future,       std::packaged_task
#include &lt;functional&gt;           // std::function,     std::bind
#include &lt;utility&gt;              // std::move,         std::forward

class ThreadPool {
public:
    inline ThreadPool(size_t threads) : stop(false) {
        for(size_t i = 0;i&lt;threads;&#43;&#43;i)
            workers.emplace_back([this] {
                for(;;) {
                    std::function&lt;void()&gt; task;
                    {
                        std::unique_lock&lt;std::mutex&gt; lock(this-&gt;queue_mutex);
                        this-&gt;condition.wait(lock,
                            [this]{ return this-&gt;stop || !this-&gt;tasks.empty(); });
                        if(this-&gt;stop &amp;&amp; this-&gt;tasks.empty())
                            return;
                        task = std::move(this-&gt;tasks.front());
                        this-&gt;tasks.pop();
                    }
                    task();
                }
            });
    }
    template&lt;class F, class... Args&gt;
    auto enqueue(F&amp;&amp; f, Args&amp;&amp;... args)
    -&gt; std::future&lt;typename std::result_of&lt;F(Args...)&gt;::type&gt; {
        using return_type = typename std::result_of&lt;F(Args...)&gt;::type;
        auto task = std::make_shared&lt; std::packaged_task&lt;return_type()&gt;&gt;(
            std::bind(std::forward&lt;F&gt;(f), std::forward&lt;Args&gt;(args)...)
        );
        std::future&lt;return_type&gt; res = task-&gt;get_future();
        {
            std::unique_lock&lt;std::mutex&gt; lock(queue_mutex);
            if(stop)
                throw std::runtime_error(&#34;enqueue on stopped ThreadPool&#34;);
            tasks.emplace([task]{ (*task)(); });
        }
        condition.notify_one();
        return res;
    }
    inline ~ThreadPool() {
        {
            std::unique_lock&lt;std::mutex&gt; lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for(std::thread &amp;worker: workers)
            worker.join();
    }
private:
    std::vector&lt;std::thread&gt; workers;
    std::queue&lt;std::function&lt;void()&gt;&gt; tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};
#endif /* ThreadPool_hpp */
```

where C&#43;&#43;11 features are:

1. lambda expression
2. std::thread
3. std::mutex, std::unique_lock
4. std::condition_variable
5. std::future, std::packaged_task
6. std::function, std::bind
7. std::shared_ptr, std::make_shared
8. std::move, std::forward

It is worth mentioning that this piece of code is very good and classic. 
It is very suitable for inclusion in our own well-preserved code snippet, and can be used immediately when it is necessary to use the thread pool for development.

## References

1. [Thread pool - wikipedia](https://en.wikipedia.org/wiki/Thread_pool)
2. [RAII - wikipedia](https://en.wikipedia.org/wiki/RAII)
3. [Monitor (Synchronization) - wikipedia](https://en.wikipedia.org/wiki/Monitor_%28synchronization%29)
4. [C&#43;&#43; Standard Library Reference](http://en.cppreference.com/w/Main_Page)
5. [C&#43;&#43; Concurrency in Action](https://www.manning.com/books/c-plus-plus-concurrency-in-action)


### Copyright

The code is modified from this project: [https://github.com/progschj/ThreadPool](https://github.com/progschj/ThreadPool).
